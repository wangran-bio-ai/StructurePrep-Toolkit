#!/usr/bin/env python3
"""contacts_agent.py — quantify interface *rearrangement* (not just contact loss).

Per-residue contact occupancy + peptide sliding distance + peptide reorientation.

Usage:
    python contacts_agent.py WT Y262A W282A E284P S261E [--ref WT]

Outputs (under ANALYSIS/):
    <V>/contacts/resid_occupancy.csv   per-residue contact occupancy
    contacts_compare.csv               per-residue difference vs reference
    contacts_compare.png               top: occupancy curves; bottom: difference
    peptide_motion.png                 peptide COM sliding + orientation change
    peptide_motion_summary.csv         one row per variant
"""
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from spconfig import cfg as CFG

B = "/root/autodl-tmp/GROMACS"
GMX = "/root/gromacs_cuda/bin/gmx"


def log(m):
    print(time.strftime("[%H:%M:%S]"), m, flush=True)


def paths(v):
    d = f"{B}/{v}/out_biop" if v == "WT" else f"{B}/{v}/out_biop_{v}"
    return d, f"{B}/ANALYSIS/{v}"


def ensure_dry(v, dt=200):
    """Protein + PEB only, whole molecules, centered: needed so contact/motion are physical."""
    d, an = paths(v)
    tpr, ndx, xtc = f"{an}/dry.tpr", f"{an}/dry.ndx", f"{an}/dry_full.xtc"
    if os.path.exists(tpr) and os.path.exists(xtc):
        return tpr, ndx, xtc
    if not os.path.exists(ndx):
        subprocess.run(f'{GMX} select -s {d}/md.tpr -select \'group "Protein" or group "PEB"\' -on {ndx}',
                       shell=True, capture_output=True, text=True)
    if not os.path.exists(tpr):
        subprocess.run(f"printf '0\\n' | {GMX} convert-tpr -s {d}/md.tpr -n {ndx} -o {tpr}",
                       shell=True, capture_output=True, text=True)
    traj = [f for f in os.listdir(d) if f.startswith("md_") and f.endswith(".xtc")][0]
    tmp = f"{an}/tmp_nojump_full.xtc"
    subprocess.run(f"printf '0\\n0\\n' | {GMX} trjconv -s {d}/md.tpr -f {d}/{traj} -n {ndx} "
                   f"-o {tmp} -pbc nojump -dt {dt}", shell=True, capture_output=True, text=True)
    subprocess.run(f"printf '0\\n0\\n' | {GMX} trjconv -s {d}/md.tpr -f {tmp} -n {ndx} "
                   f"-o {xtc} -pbc mol -center -ur compact", shell=True, capture_output=True, text=True)
    if os.path.exists(tmp):
        os.remove(tmp)
    return tpr, ndx, xtc


def residue_map(v):
    d, _ = paths(v)
    L = open(f"{d}/complex.gro").read().splitlines()
    n = int(L[1])
    res, cur = [], None
    for ln in L[2:2 + n]:
        rid, rname = int(ln[0:5]), ln[5:10].strip()
        if cur is None or cur[0] != rid or cur[1] != rname:
            cur = (rid, rname, 0)
            res.append(cur)
        res[-1] = (cur[0], cur[1], res[-1][2] + 1)
    return [r for r in res if r[1] != "PEB"]


def coords(v):
    an = paths(v)[1]
    out = f"{an}/dry_traj.xvg"
    tpr, ndx, xtc = ensure_dry(v)
    if not os.path.exists(out):
        log(f"{v}: dumping coordinates")
        subprocess.run(f"printf '0\\n' | {GMX} traj -s {tpr} -f {xtc} -n {ndx} -ox {out}",
                       shell=True, capture_output=True, text=True)
    rows = [l.split() for l in open(out) if l.strip() and l[0] not in "#@"]
    arr = np.array(rows, dtype=float)
    return arr[:, 0] / 1000.0, arr[:, 1:].reshape(len(arr), -1, 3)


def analyze(v, cut=None):
    cut = CFG.getfloat("contacts", "cutoff_nm", 0.4) if cut is None else cut
    prot = residue_map(v)
    t, X = coords(v)
    nprot = sum(r[2] for r in prot)
    P, L = X[:, :nprot, :], X[:, nprot:, :]
    starts = np.cumsum([0] + [r[2] for r in prot])

    # per-residue contact occupancy
    occ = np.zeros(len(prot))
    for i, (rid, rname, nat) in enumerate(prot):
        s = starts[i]
        d = np.linalg.norm(P[:, s:s + nat, None, :] - L[:, None, :, :], axis=-1).min(axis=(1, 2))
        occ[i] = np.mean(d < cut)
    df = pd.DataFrame({"resid": [r[0] for r in prot], "resname": [r[1] for r in prot],
                       "occupancy": occ})
    an = paths(v)[1]
    os.makedirs(f"{an}/contacts", exist_ok=True)
    df.to_csv(f"{an}/contacts/resid_occupancy.csv", index=False)

    # peptide position RELATIVE to the protein, in the protein's own frame:
    # first best-fit superimpose the protein (CA atoms) onto frame 0, then take the
    # peptide COM. Without this, rigid-body rotation of the whole complex would be
    # mistaken for peptide sliding.
    ca_idx = [starts[i] + 1 for i in range(len(prot))]
    ca = P[:, ca_idx, :]
    ref_ca = ca[0] - ca[0].mean(axis=0)
    rel = np.zeros((len(t), 3))
    axes = np.zeros((len(t), 3))
    pep_vec0 = L[0, 28, :] - L[0, 0, :]
    for k in range(len(t)):
        mov = ca[k] - ca[k].mean(axis=0)
        H = mov.T @ ref_ca
        U, S, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        R = Vt.T @ np.diag([1, 1, d]) @ U.T
        rel[k] = R @ (L[k].mean(axis=0) - ca[k].mean(axis=0))
        axes[k] = R @ (L[k, 28, :] - L[k, 0, :])
    disp = np.linalg.norm(rel - rel[0], axis=1)
    axes = axes / np.linalg.norm(axes, axis=1, keepdims=True)

    idx = {r[0]: i for i, r in enumerate(prot)}
    if 262 in idx and 286 in idx:
        v262 = P[:, starts[idx[262]] + 1, :]
        v286 = P[:, starts[idx[286]] + 1, :]
        ref = v286 - v262
        ref = ref / np.linalg.norm(ref, axis=1, keepdims=True)
        ang = np.degrees(np.arccos(np.clip(np.sum(axes * ref, axis=1), -1, 1)))
        ang = ang - ang[0]
    else:
        ang = np.zeros(len(t))

    motion = pd.DataFrame({"t_ns": t, "com_disp_nm": disp, "axis_angle_deg": ang})
    motion.to_csv(f"{an}/contacts/peptide_motion.csv", index=False)
    log(f"{v}: COM sliding {disp.mean():.2f} nm (max {disp.max():.2f}) | "
        f"|angle| {np.abs(ang).mean():.1f} deg | frames {len(t)}")
    return df, motion


def main():
    args = sys.argv[1:]
    ref = "WT"
    if "--ref" in args:
        i = args.index("--ref")
        ref = args[i + 1]
        del args[i:i + 2]
    variants = args or ["WT", "Y262A"]
    if ref not in variants:
        variants = [ref] + variants

    occ, mot = {}, {}
    for v in variants:
        occ[v], mot[v] = analyze(v)

    rows = []
    for v in variants:
        rows.append({
            "variant": v,
            "sum_occupancy": occ[v].occupancy.sum(),
            "n_core_residues": int((occ[v].occupancy > CFG.getfloat("contacts", "core_occupancy", 0.5)).sum()),
            "mean_com_disp_nm": mot[v].com_disp_nm.mean(),
            "mean_abs_angle_deg": float(np.abs(mot[v].axis_angle_deg).mean()),
            "d_sum_occupancy_vs_ref": occ[v].occupancy.sum() - occ[ref].occupancy.sum(),
        })
    summ = pd.DataFrame(rows)
    summ.to_csv(f"{B}/ANALYSIS/peptide_motion_summary.csv", index=False)
    print("\n=== interface rearrangement summary ===")
    print(summ.to_string(index=False))

    if ref in occ:
        cmp = pd.concat([
            pd.DataFrame({"variant": v,
                          "resid": occ[v].resid.values,
                          "d_occ": occ[v].occupancy.values - occ[ref].occupancy.values})
            for v in variants if v != ref])
        cmp.to_csv(f"{B}/ANALYSIS/contacts_compare.csv", index=False)

        fig, ax = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
        ax[0].plot(occ[ref].resid, occ[ref].occupancy, lw=1.2, label=ref, color="k")
        for v in variants:
            if v != ref:
                ax[0].plot(occ[v].resid, occ[v].occupancy, lw=.9, alpha=.75, label=v)
        ax[0].set(ylabel="contact occupancy", title="per-residue contact occupancy (protein <-> PEB, <4 A)")
        ax[0].legend(fontsize=8)
        for v in variants:
            if v == ref:
                continue
            d = cmp[cmp.variant == v]
            ax[1].bar(d.resid, d.d_occ, width=1.0, alpha=.6, label=f"{v} - {ref}")
        ax[1].axhline(0, c="k", lw=.8)
        ax[1].set(xlabel="residue", ylabel=f"delta occupancy vs {ref}", title="interface rearrangement")
        ax[1].legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(f"{B}/ANALYSIS/contacts_compare.png", dpi=150)
        plt.close()

        fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
        for v in variants:
            ax[0].plot(mot[v].t_ns, mot[v].com_disp_nm, lw=.9, label=v)
            ax[1].plot(mot[v].t_ns, mot[v].axis_angle_deg, lw=.9, label=v)
        ax[0].set(xlabel="time (ns)", ylabel="peptide COM displacement (nm)", title="peptide sliding")
        ax[1].set(xlabel="time (ns)", ylabel="axis angle change (deg)", title="peptide reorientation")
        for a in ax:
            a.legend(fontsize=7)
        plt.tight_layout()
        plt.savefig(f"{B}/ANALYSIS/peptide_motion.png", dpi=150)
        plt.close()
        log("written: contacts_compare.png / peptide_motion.png / peptide_motion_summary.csv")


if __name__ == "__main__":
    main()
