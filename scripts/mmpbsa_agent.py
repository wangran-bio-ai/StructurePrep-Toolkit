#!/usr/bin/env python3
"""mmpbsa_agent.py — MM/GBSA（或 MM/PBSA）批量分析

用法:
  python mmpbsa_agent.py WT Y262A S261E --ref WT            # 从轨迹重算
  python mmpbsa_agent.py WT Y262A S261E --ref WT --reuse    # 已有 .dat 就只汇总
  python mmpbsa_agent.py WT --window 50 100 --frames 125

流程（每一步都可跳过，结果缓存在 ANALYSIS/<V>/mmpbsa/）:
  1. 只保留 蛋白+配体 -> dry.tpr（convert-tpr）
  2. 先 -pbc nojump 再 -pbc mol -center（**不做 nojump 会让配体被成像到盒子另一侧，ΔG≈0**）
  3. gmx_MMPBSA GB(igb=5) 或 PB
  4. 解析 .dat -> 均值/SD/SEM、ΔΔG vs 参考、CSV + 条形图

判读边界（写论文时照抄）:
  * ΔEEL 与 ΔEGB/ΔEPB 必须成对出现，单独引用 ΔEEL 会得出错误结论
  * 单轨迹近似、无熵项 => 只报 ΔΔG 的方向与量级，不报绝对 ΔG
  * 若两个变体在 MD 里都保持结合，"结合姿态下的 ΔG 相近"是正常结果，不代表亲和力相同
"""
import argparse
import os
import re
import subprocess
import sys
import time

import numpy as np

B = "/root/autodl-tmp/GROMACS"
GMX = "/root/gromacs_cuda/bin/gmx"
AMBER_ENV = "/root/autodl-tmp/mmpbsa"
COLS = ["VDWAALS", "EEL", "EGB", "ESURF", "TOTAL"]


def log(m):
    print(time.strftime("[%H:%M:%S]"), m, flush=True)


def sh(c):
    r = subprocess.run(["bash", "-lc", c], capture_output=True, text=True)
    return (r.stdout + r.stderr).strip()


def paths(v):
    d = f"{B}/{v}/out_biop" if v == "WT" else f"{B}/{v}/out_biop_{v}"
    return d, f"{B}/ANALYSIS/{v}"


def prep(v, t0, t1, dt_ps):
    d, an = paths(v)
    mdir = f"{an}/mmpbsa"
    os.makedirs(mdir, exist_ok=True)
    ndx, tpr, xtc = f"{an}/dry.ndx", f"{an}/dry.tpr", f"{an}/dry_{int(t0)}-{int(t1)}ns.xtc"
    gndx = f"{an}/dry_groups.ndx"
    if not os.path.exists(ndx):
        sh(f'{GMX} select -s {d}/md.tpr -select \'group "Protein" or group "PEB"\' -on {ndx}')
    if not os.path.exists(tpr):
        sh(f"printf '0\\n' | {GMX} convert-tpr -s {d}/md.tpr -n {ndx} -o {tpr}")
    if not os.path.exists(xtc):
        traj = [f for f in os.listdir(d) if f.startswith("md_") and f.endswith(".xtc")][0]
        tmp = f"{an}/tmp_nojump.xtc"
        sh(f"printf '0\\n0\\n' | {GMX} trjconv -s {d}/md.tpr -f {d}/{traj} -n {ndx} "
           f"-o {tmp} -pbc nojump -b {int(t0*1000)} -e {int(t1*1000)} -dt {dt_ps}")
        sh(f"printf '0\\n0\\n' | {GMX} trjconv -s {d}/md.tpr -f {tmp} -n {ndx} "
           f"-o {xtc} -pbc mol -center -ur compact")
        if os.path.exists(tmp):
            os.remove(tmp)
    if not os.path.exists(gndx):
        sh(f"cd {an} && printf 'q\\n' | {GMX} make_ndx -f dry.tpr -o dry_groups.ndx")
    g = {}
    for ln in open(gndx):
        s = ln.strip()
        if s.startswith("["):
            g[s.strip("[] \t")] = len(g)
    return d, an, tpr, gndx, xtc, g


def run_one(v, t0, t1, dt_ps, nframes, method, salt, igb):
    d, an, tpr, ndx, xtc, g = prep(v, t0, t1, dt_ps)
    mdir = f"{an}/mmpbsa"
    total_frames = max(2, int((t1 - t0) * 1000 / dt_ps))
    interval = max(1, total_frames // nframes)
    inp = f"{mdir}/mmpbsa.in"
    with open(inp, "w") as f:
        f.write(f'&general\n  sys_name="{v}", startframe=1, endframe={total_frames}, '
                f'interval={interval}, verbose=2,\n/\n')
        if method == "gb":
            f.write(f'&gb\n  igb={igb}, saltcon={salt},\n/\n')
        else:
            f.write(f'&pb\n  istrng={salt}, fillratio=4.0,\n/\n')
    env = dict(os.environ)
    env["PATH"] = f"{AMBER_ENV}/bin:{os.path.dirname(GMX)}:" + env.get("PATH", "")
    env["AMBERHOME"] = AMBER_ENV
    cmd = [f"{AMBER_ENV}/bin/gmx_MMPBSA", "-O", "-i", inp,
           "-cs", tpr, "-ci", ndx, "-cg", str(g["Protein"]), str(g["PEB"]),
           "-ct", xtc, "-cp", f"{d}/topol.top",
           "-o", f"{mdir}/{v}.csv", "-eo", f"{mdir}/{v}.dat", "-prefix", f"{v}_"]
    logf = f"{mdir}/{v}.log"
    with open(logf, "w") as fh:
        rc = subprocess.run(cmd, cwd=mdir, stdout=fh, stderr=subprocess.STDOUT, env=env).returncode
    # gmx_MMPBSA 末尾调用 PyQt5 GUI 会返回 1，但计算已完成 —— 以 .dat 是否存在为准
    ok = os.path.exists(f"{mdir}/{v}.dat")
    log(f"{v}: rc={rc} dat={ok}" + ("" if ok else "  tail=" + sh(f"tail -n 5 {logf}")[-300:]))
    return ok


def parse(v, method):
    _, an = paths(v)
    p = None
    for cand in (f"{an}/mmpbsa/{v}.dat", f"{an}/mmpbsa/{v}_full.dat"):
        if os.path.exists(cand):
            p = cand
            break
    if p is None:
        return None
    blk = open(p, errors="ignore").read().split("Delta Energy Terms")[-1]
    hdr = [l for l in blk.splitlines() if l.startswith("Frame #")][0].split(",")
    cols = [c for c in COLS if c in hdr] + (["EPB", "TOTAL"] if "EPB" in hdr else [])
    rows = [l for l in blk.splitlines() if re.match(r"^\s*\d+,", l)]
    arr = np.array([[float(x) for x in l.split(",")] for l in rows])
    idx = {c: hdr.index(c) for c in hdr if c != "Frame #"}
    return {c: (arr[:, idx[c]].mean(), arr[:, idx[c]].std(), len(arr)) for c in idx}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("variants", nargs="*", default=["WT"])
    ap.add_argument("--ref", default="WT")
    ap.add_argument("--window", nargs=2, type=float, default=[50, 100])
    ap.add_argument("--frames", type=int, default=125)
    ap.add_argument("--dt", type=float, default=200.0, help="抽帧间隔 ps")
    ap.add_argument("--method", choices=["gb", "pb"], default="gb")
    ap.add_argument("--salt", type=float, default=0.15)
    ap.add_argument("--igb", type=int, default=5)
    ap.add_argument("--reuse", action="store_true")
    a = ap.parse_args()

    variants = a.variants or ["WT"]
    for v in variants:
        if not a.reuse:
            run_one(v, a.window[0], a.window[1], a.dt, a.frames, a.method, a.salt, a.igb)

    res = {}
    for v in variants:
        st = parse(v, a.method)
        if st:
            res[v] = st
    if not res:
        sys.exit("no results")

    key = "TOTAL" if "TOTAL" in next(iter(res.values())) else "TOTAL "
    print("\n=== ΔG_bind（%s, kcal/mol, n 帧）===" % a.method.upper())
    for v, st in res.items():
        m, s, n = st["TOTAL"]
        sem = s / np.sqrt(n)
        print(f"{v:8s} ΔG = {m:7.2f} ± {sem:4.2f} (SEM, n={n})   "
              f"[VDW {st['VDWAALS'][0]:.1f} | EEL {st['EEL'][0]:.1f} | "
              f"{(st.get('EGB') or st.get('EPB'))[0]:.1f} | ESURF {st['ESURF'][0]:.1f}]")

    ref = a.ref if a.ref in res else next(iter(res))
    print(f"\n=== ΔΔG vs {ref}（正值 = 结合变弱）===")
    rows = []
    for v, st in res.items():
        m, s, n = st["TOTAL"]
        rm, rs, rn = res[ref]["TOTAL"]
        d = m - rm
        err = (s ** 2 / n + rs ** 2 / rn) ** 0.5
        rows.append((v, m, s / np.sqrt(n), d, err))
        if v != ref:
            print(f"{v:8s} ΔΔG = {d:+6.2f} ± {err:.2f} kcal/mol")

    an = f"{B}/ANALYSIS"
    out = f"{an}/mmpbsa_summary.csv"
    with open(out, "w") as f:
        f.write("variant,dG_kcal,sem,n_frames,dDG_vs_%s,dDG_err\n" % ref)
        for v, m, sem, d, err in rows:
            f.write(f"{v},{m:.3f},{sem:.3f},{res[v]['TOTAL'][2]},{d:.3f},{err:.3f}\n")
    log(f"summary -> {out}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        vs = [r[0] for r in rows]
        gs = [r[1] for r in rows]
        es = [r[2] for r in rows]
        fig, ax = plt.subplots(figsize=(5.2, 3.4))
        ax.bar(vs, gs, yerr=es, capsize=4, color="C0")
        ax.set(ylabel="ΔG_bind (kcal/mol)", title=f"MM-{a.method.upper()} (single trajectory)")
        ax.axhline(0, c="k", lw=.8)
        plt.tight_layout()
        figp = f"{an}/mmpbsa_bar.png"
        plt.savefig(figp, dpi=150)
        log(f"figure -> {figp}")
    except Exception as e:
        log(f"plot skipped: {e}")


if __name__ == "__main__":
    main()
