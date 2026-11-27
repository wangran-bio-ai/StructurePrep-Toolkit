#!/usr/bin/env python3
"""amber_decomp.py — 用 AmberTools 的 MMPBSA.py 做逐残基分解（gmx_MMPBSA 的 decomp 解析器在本体系失效）

流程：
  1. dry 轨迹(50–100 ns，每 2 ns 一帧) -> 多模型 PDB -> cpptraj 转 mdcrd
  2. MMPBSA.py -cp COM.prmtop -rp REC.prmtop -lp LIG.prmtop -y traj.mdcrd
  3. 解析 FINAL_DECOMP_MMPBSA.dat 的每残基 DELTA TOTAL
用法: python amber_decomp.py WT W282A R279E S261E
"""
import os
import re
import subprocess
import sys

import pandas as pd

B = "/root/autodl-tmp/GROMACS"
ENV = "/root/autodl-tmp/mmpbsa"
GMX = "/root/gromacs_cuda/bin/gmx"
OUT = f"{B}/ANALYSIS"


def sh(c, tag=""):
    r = subprocess.run(["bash", "-c", c], capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[{tag}] rc={r.returncode}: {(r.stdout + r.stderr)[-300:]}")
    return r.stdout + r.stderr


def prep(v, frames_step_ps=2000):
    an = f"{B}/ANALYSIS/{v}"
    mdir = f"{an}/mmpbsa"
    pdb = f"{mdir}/traj_frames.pdb"
    mdcrd = f"{mdir}/traj.mdcrd"
    if os.path.exists(mdcrd) and os.path.getsize(mdcrd) > 10000:
        return mdir, mdcrd
    subprocess.run([GMX, "trjconv", "-s", f"{an}/dry.tpr", "-f", f"{an}/dry_50-100ns.xtc",
                    "-n", f"{an}/dry_groups.ndx", "-o", pdb, "-dt", str(frames_step_ps)],
                   input="0\n", capture_output=True, text=True)
    env = dict(os.environ, PATH=f"{ENV}/bin:" + os.environ["PATH"], AMBERHOME=ENV)
    subprocess.run([f"{ENV}/bin/cpptraj", "-p", f"{mdir}/COM.prmtop", "-y", pdb,
                    "-x", mdcrd], capture_output=True, text=True, env=env)
    print(f"[{v}] frames -> {os.path.basename(mdcrd)} "
          f"({os.path.getsize(mdcrd)//1024 if os.path.exists(mdcrd) else 0} KB)")
    return mdir, mdcrd


def run(v):
    mdir, mdcrd = prep(v)
    inp = f"{mdir}/amber_decomp.in"
    with open(inp, "w") as f:
        f.write("&general\n  startframe=1, endframe=999999, interval=1,\n/\n")
        f.write("&gb\n  igb=5, saltcon=0.15,\n/\n")
        f.write("&decomp\n  idecomp=1, print_res=\"all\", dec_verbose=1,\n/\n")
    env = dict(os.environ, PATH=f"{ENV}/bin:" + os.environ["PATH"], AMBERHOME=ENV)
    cmd = [f"{ENV}/bin/MMPBSA.py", "-O", "-i", inp,
           "-o", f"{mdir}/amber_{v}.dat",
           "-cp", f"{mdir}/COM.prmtop", "-rp", f"{mdir}/REC.prmtop",
           "-lp", f"{mdir}/LIG.prmtop", "-y", mdcrd,
           "-eo", f"{mdir}/amber_{v}_perframe.csv",
           "-do", f"{mdir}/amber_{v}_decomp.csv"]
    print(f"[{v}] running MMPBSA.py ...")
    r = subprocess.run(cmd, cwd=mdir, capture_output=True, text=True, env=env)
    ok = os.path.exists(f"{mdir}/amber_{v}_decomp.csv")
    print(f"[{v}] rc={r.returncode} decomp_csv={ok}")
    if not ok:
        print((r.stdout + r.stderr)[-600:])
    return f"{mdir}/amber_{v}_decomp.csv"


def parse(path):
    """Amber 的 decomp CSV：列 = Residue, Location, Internal(3), vdW(3), Elec(3),
    Polar(3), NonPolar(3), TOTAL(3)；TOTAL 平均值在第 17 列（0-based）。
    残基号是 Amber 内部重编号 1..126（对应蛋白 204..329）。"""
    out = {}
    for ln in open(path, errors="ignore"):
        parts = [p.strip() for p in ln.split(",")]
        if len(parts) < 18:
            continue
        m = re.match(r"^([A-Z]{3})\s+(\d+)$", parts[0])
        if not m:
            continue
        try:
            out[int(m.group(2)) + 203] = dict(
                TOTAL=float(parts[17]), VDW=float(parts[5]), ELE=float(parts[8]),
                GB=float(parts[11]), SA=float(parts[14]))
        except (ValueError, IndexError):
            continue
    return out or None


def main():
    variants = sys.argv[1:] or ["WT", "W282A"]
    data = {}
    for v in variants:
        p = run(v)
        if os.path.exists(p):
            d = parse(p)
            if d:
                data[v] = d
                print(f"[{v}] 残基数 {len(d)}")
    if len(data) >= 2:
        ref = variants[0]
        for v in variants[1:]:
            if v not in data:
                continue
            keys = sorted(set(data[ref]) & set(data[v]))
            df = pd.DataFrame({"resid": keys,
                               "WT_total": [data[ref][k]["TOTAL"] for k in keys],
                               f"{v}_total": [data[v][k]["TOTAL"] for k in keys],
                               "WT_vdw": [data[ref][k]["VDW"] for k in keys],
                               f"{v}_vdw": [data[v][k]["VDW"] for k in keys],
                               "WT_ele": [data[ref][k]["ELE"] for k in keys],
                               f"{v}_ele": [data[v][k]["ELE"] for k in keys]})
            df["d_total"] = df[f"{v}_total"] - df["WT_total"]
            df["d_vdw"] = df[f"{v}_vdw"] - df["WT_vdw"]
            df["d_ele"] = df[f"{v}_ele"] - df["WT_ele"]
            df["delta"] = df["d_total"]
            df = df.sort_values("delta")
            df.to_csv(f"{OUT}/decomp_{v}_vs_{ref}.csv", index=False)
            print(f"\n=== {v} vs {ref} 最负（结合更强）===")
            print(df.head(5).to_string(index=False))
            print(f"=== {v} vs {ref} 最正（结合更弱）===")
            print(df.tail(5).to_string(index=False))
            print("合计:", round(df.delta.sum(), 2))


if __name__ == "__main__":
    main()
