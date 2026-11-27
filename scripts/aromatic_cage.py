#!/usr/bin/env python3
"""aromatic_cage.py — 芳香笼 / "盖子" 的几何量化

用已有 dry 轨迹（ANALYSIS/<V>/dry_traj.xvg，501 帧）计算：
  1. Trp282 吲哚环与生物素噻吩环的**质心距离**（盖子闭合程度）
  2. 两环的**二面角**（0°=平行堆积，90°=edge-to-face）
  3. Trp282 与 Tyr262 / His259 / Phe285 的质心距离（笼子尺寸）
  4. Trp282 环平面法向与"环心→配体质心"方向的夹角（盖子倾角）

输出：ANALYSIS/cage_geometry.csv、cage_geometry.png
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

B = "/root/autodl-tmp/GROMACS"
OUT = f"{B}/ANALYSIS"
VARS = ["WT", "Y262A", "W282A", "E284P", "S261E", "R279E"]

RINGS = {
    "W282": (282, ["CG", "CD1", "CD2", "CE2", "NE1", "CE3", "CZ2", "CZ3", "CH2"]),
    "Y262": (262, ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]),
    "H259": (259, ["CG", "ND1", "CD2", "CE1", "NE2"]),
    "F285": (285, ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]),
}


def read_gro(v):
    d = f"{B}/{v}/out_biop" if v == "WT" else f"{B}/{v}/out_biop_{v}"
    L = open(f"{d}/complex.gro").read().splitlines()
    n = int(L[1])
    atoms = []
    for ln in L[2:2 + n]:
        atoms.append(dict(resid=int(ln[0:5]), resname=ln[5:10].strip(),
                          name=ln[10:15].strip()))
    return atoms


def load_traj(v):
    p = f"{OUT}/{v}/dry_traj.xvg"
    rows = [l.split() for l in open(p) if l.strip() and l[0] not in "#@"]
    arr = np.array(rows, dtype=float)
    return arr[:, 0] / 1000.0, arr[:, 1:].reshape(len(arr), -1, 3)


def ring_indices(atoms, resid, names):
    idx = [i for i, a in enumerate(atoms) if a["resid"] == resid and a["name"] in names]
    return idx


def biotin_ring_indices(atoms, X0):
    """生物素噻吩环：S 原子 + 与 S 成键(<=1.9 Å)的碳"""
    s_idx = [i for i, a in enumerate(atoms) if a["resname"] == "PEB" and a["name"].strip() == "S"]
    if not s_idx:
        s_idx = [i for i, a in enumerate(atoms)
                 if a["resname"] == "PEB" and a["name"].strip().upper().startswith("S")]
    if not s_idx:
        return None
    s = s_idx[0]
    d = np.linalg.norm(X0 - X0[s], axis=1)
    # 只取同一分子（PEB）内的原子，否则会把界面上相邻的蛋白原子也算进环里
    ring = [i for i in np.where(d < 2.0)[0]
            if i != s and atoms[i]["resname"] == "PEB"]
    return [s] + ring


def centroid_plane(X, idx):
    P = X[:, idx, :].mean(axis=1)                 # (nframe, 3)
    Q = X[:, idx, :] - P[:, None, :]
    normals = []
    for k in range(X.shape[0]):
        _, _, vt = np.linalg.svd(Q[k], full_matrices=False)
        normals.append(vt[2])
    return P, np.array(normals)


def angle_between(n1, n2):
    c = np.abs(np.sum(n1 * n2, axis=1)) / (np.linalg.norm(n1, axis=1) * np.linalg.norm(n2, axis=1))
    return np.degrees(np.arccos(np.clip(c, 0, 1)))


def main():
    rows = []
    for v in VARS:
        atoms = read_gro(v)
        t, X = load_traj(v)
        X0 = X[0]
        biotin = biotin_ring_indices(atoms, X0)
        if biotin is None:
            print(f"[{v}] 找不到生物素 S 原子"); continue

        cP, cN = centroid_plane(X, biotin)
        rings = {}
        for tag, (resid, names) in RINGS.items():
            idx = ring_indices(atoms, resid, names)
            if len(idx) < 3:
                continue
            rings[tag] = centroid_plane(X, idx)

        if "W282" in rings:
            wP, wN = rings["W282"]
        else:
            ca = [i for i, a in enumerate(atoms) if a["resid"] == 282 and a["name"] == "CA"]
            wP = X[:, ca, :].mean(axis=1)
            wN = np.tile(np.array([0.0, 0.0, 1.0]), (X.shape[0], 1))
        d_lig = np.linalg.norm(cP - wP, axis=1) * 10.0   # nm -> Å
        ang_lig = angle_between(cN, wN) if "W282" in rings else np.full(len(d_lig), np.nan)

        # 盖子倾角：环法向 vs 环心→配体质心
        vec = cP - wP
        cos = np.abs(np.sum(vec * wN, axis=1)) / (np.linalg.norm(vec, axis=1) + 1e-9)
        tilt = np.degrees(np.arccos(np.clip(cos, 0, 1)))

        rec = dict(variant=v, n_frames=len(t),
                   W282_biotin_dist_A=round(d_lig.mean(), 2),
                   W282_biotin_dist_sd=round(d_lig.std(), 2),
                   W282_biotin_angle_deg=round(ang_lig.mean(), 1),
                   lid_tilt_deg=(round(tilt.mean(), 1) if "W282" in rings else None))
        for tag in ("Y262", "H259", "F285"):
            if tag in rings:
                d = np.linalg.norm(rings[tag][0] - wP, axis=1) * 10.0
                rec[f"W282_{tag}_dist_A"] = round(d.mean(), 2)
        rows.append(rec)
        print(f"[{v}] W282–biotin {d_lig.mean():.2f} Å, 环二面角 {ang_lig.mean():.0f}°, "
              f"盖子倾角 {tilt.mean():.0f}°")

    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/cage_geometry.csv", index=False)
    print(df.to_string(index=False))

    # 图：距离与角度分布
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.8))
    sub = df.set_index("variant")
    axes[0].bar(sub.index, sub.W282_biotin_dist_A, yerr=sub.W282_biotin_dist_sd, capsize=4)
    axes[0].set(ylabel="W282–biotin centroid (Å)", title="aromatic lid closure")
    axes[1].bar(sub.index, sub.W282_biotin_angle_deg, color="C1")
    axes[1].set(ylabel="interplanar angle (deg)", title="ring stacking geometry")
    axes[2].bar(sub.index, sub.lid_tilt_deg, color="C2")
    axes[2].set(ylabel="lid tilt (deg)", title="Trp282 ring orientation")
    for a in axes:
        a.tick_params(axis="x", labelrotation=45)
    plt.tight_layout()
    plt.savefig(f"{OUT}/cage_geometry.png", dpi=150)
    print("图 ->", f"{OUT}/cage_geometry.png")


if __name__ == "__main__":
    main()
