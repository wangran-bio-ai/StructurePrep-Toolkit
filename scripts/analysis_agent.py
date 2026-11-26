#!/usr/bin/env python3
"""MD 界面分析智能体 (analysis_agent.py)

一次跑完所有变体，输出：
  ANALYSIS/<V>/{xvg,logs,figs,summary.csv}
  ANALYSIS/comparison/{metrics.csv, ranking.csv, figs/*.png}

用法:
  python analysis_agent.py                 # 默认分析 WT Y262A W282A（存在的才跑）
  python analysis_agent.py WT R279E --ref WT
  python analysis_agent.py --reuse         # 已有 xvg 就复用，只做汇总

设计要点（来自 Y262A 的实战教训）：
  1. 目录命名不统一：WT 用 out_biop，其它用 out_biop_<V> → 自动探测
  2. gmx 路径不固定 → 自动探测候选路径
  3. PBC 跳变会让肽位移出现 8-12 nm 假尖峰 → 用中位数 + 剔除离群帧
  4. PEB 是 GAFF 原子名，tpr 无元素类型 → gmx hbond 不可用，改用几何判据
  5. 接触占时比只能当"闸门"，不能反推结合强弱
"""
import os
import sys
import glob
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/root/autodl-tmp/GROMACS"
AN = f"{BASE}/ANALYSIS"
GMX_CANDS = ["/root/gromacs_cuda/bin/gmx",
             "/root/gromacs_cuda_244/bin/gmx",
             "/root/autodl-tmp/md/bin/gmx"]
PBC_CUT = 2.0          # nm，肽位移超过它视为 PBC 跳变
TRIM_NS = 50.0         # 默认分析窗口终点（ns）：WT 100 ns、突变体 50 ns，统一截到 0–50 ns
TRIM_FROM = 0.0        # 窗口起点（ns）。>0 时输出到带后缀的目录，且不参与跨变体对比
TRIM_TO = TRIM_NS
DEFAULT_VARS = ["WT", "Y262A", "W282A"]
REF = "WT"


def tw():
    """轨迹窗口。用 --from/--to 覆盖时为"窗口一致性检验"，不与其他变体混比。"""
    return ["-b", f"{TRIM_FROM * 1000:.0f}", "-e", f"{TRIM_TO * 1000:.0f}"]


# ---------------------------------------------------------------- 基础工具
def find_gmx():
    for c in GMX_CANDS:
        if os.path.exists(c):
            return c
    hit = glob.glob("/root/*/gromacs_cuda/bin/gmx") + glob.glob("/root/*/bin/gmx")
    if hit:
        return hit[0]
    raise SystemExit("gmx not found")


GMX = find_gmx()


def var_dir(v):
    for cand in (f"{BASE}/{v}/out_biop_{v}", f"{BASE}/{v}/out_biop"):
        if os.path.isdir(cand):
            return cand
    return None


def run(args, stdin_text=None, log=None):
    r = subprocess.run([GMX] + args, input=stdin_text, text=True, capture_output=True)
    if log:
        with open(log, "w") as f:
            f.write(r.stdout + r.stderr)
    if r.returncode != 0:
        raise RuntimeError(f"gmx {args[0]} failed: " + (r.stdout + r.stderr)[-600:])
    return r


def read_xvg(path):
    rows = []
    for ln in open(path):
        s = ln.strip()
        if not s or s[0] in "#@":
            continue
        rows.append([float(x) for x in s.split()])
    return np.array(rows)


def names_from_itp(path):
    out, sec = [], None
    for ln in open(path, errors="ignore"):
        s = ln.strip()
        if s.startswith("["):
            sec = s.strip("[] \t")
            continue
        if sec == "atoms" and s and not s.startswith(";"):
            f = s.split()
            if len(f) >= 6 and f[0].isdigit():
                out.append(f[4])
    return out


# ---------------------------------------------------------------- 单变体分析
def analyze(v, reuse=False):
    D = var_dir(v)
    if D is None:
        print(f"[skip] {v}: 目录不存在")
        return None
    OUT = f"{AN}/{v}" + ("" if TRIM_FROM == 0 else f"_w{int(TRIM_FROM)}_{int(TRIM_TO)}")
    for sub in ("xvg", "logs", "figs"):
        os.makedirs(f"{OUT}/{sub}", exist_ok=True)

    TPR = f"{D}/md.tpr"
    xtcs = sorted(glob.glob(f"{D}/md_*.part0*.xtc")) or sorted(glob.glob(f"{D}/md_*.xtc"))
    XTC = xtcs[0]
    IDX = f"{OUT}/index.ndx"
    print(f"[run] {v}  dir={D}")

    # index
    if not (reuse and os.path.exists(IDX)):
        run(["make_ndx", "-f", TPR, "-o", IDX], stdin_text="q\n", log=f"{OUT}/logs/makenndx.log")
    G = {}
    for ln in open(f"{OUT}/logs/makenndx.log"):
        p = ln.split()
        if len(p) >= 4 and p[0].isdigit() and p[1] != ":":
            G[ln.split(":")[0].split(None, 1)[1].strip()] = int(p[0])
    n_prot = int([l for l in open(f"{OUT}/logs/makenndx.log") if "Protein " in l][0]
                 .split(":")[1].split()[0])

    def need(name):
        return not (reuse and os.path.exists(f"{OUT}/xvg/{name}.xvg"))

    # 1) RMSD / Rg / RMSF
    if need("rmsd_bb"):
        run(["rms", "-s", TPR, "-f", XTC, "-n", IDX, "-o", f"{OUT}/xvg/rmsd_bb.xvg", "-tu", "ns"] + tw(),
            stdin_text=f"{G['Backbone']}\n{G['Backbone']}\n", log=f"{OUT}/logs/rmsd.log")
    if need("gyrate"):
        run(["gyrate", "-s", TPR, "-f", XTC, "-n", IDX, "-o", f"{OUT}/xvg/gyrate.xvg"] + tw(),
            stdin_text=f"{G['Protein']}\n", log=f"{OUT}/logs/gyrate.log")
    if need("rmsf"):
        run(["rmsf", "-s", TPR, "-f", XTC, "-n", IDX, "-o", f"{OUT}/xvg/rmsf.xvg", "-res"] + tw(),
            stdin_text=f"{G['C-alpha']}\n", log=f"{OUT}/logs/rmsf.log")

    # 2) 接触 + 最小距离 + 肽位移
    if need("ncontacts"):
        run(["mindist", "-s", TPR, "-f", XTC, "-n", IDX, "-d", "0.4",
             "-od", f"{OUT}/xvg/mindist_peb.xvg", "-on", f"{OUT}/xvg/ncontacts.xvg"] + tw(),
            stdin_text=f"{G['Protein']}\n{G['PEB']}\n", log=f"{OUT}/logs/mindist.log")
    if need("rmsd_peb"):
        run(["rms", "-s", TPR, "-f", XTC, "-n", IDX, "-o", f"{OUT}/xvg/rmsd_peb.xvg", "-tu", "ns"] + tw(),
            stdin_text=f"{G['Protein']}\n{G['PEB']}\n", log=f"{OUT}/logs/rmsd_pep.log")

    # 3) 279 位距离分布（最近距离；WT 与突变体可比）
    if need("dist_279_peb"):
        run(["pairdist", "-s", TPR, "-f", XTC, "-ref", "resnr 279", "-sel", "resname PEB",
             "-o", f"{OUT}/xvg/dist_279_peb.xvg"] + tw(), log=f"{OUT}/logs/pairdist279.log")

    # 4) buried SASA（复合物 / 蛋白 / 肽 三次）
    if need("sasa_complex"):
        for sel, tag in ((f'group "Protein" or group "PEB"', "sasa_complex"),
                         (f'group "Protein"', "sasa_protein"),
                         (f'group "PEB"', "sasa_peb")):
            run(["sasa", "-s", TPR, "-f", XTC, "-n", IDX, "-surface", sel, "-dt", "0.5",
                 "-o", f"{OUT}/xvg/{tag}.xvg"] + tw(), log=f"{OUT}/logs/{tag}.log")

    # 5) 候选氢键（几何判据，绕开 GAFF 原子名无元素类型的问题）
    if need("polar_pairs"):
        pn = names_from_itp(f"{D}/topol.top")
        qn = names_from_itp(f"{D}/PEPB.acpype/PEPB_GMX.itp")
        prot_idx = [i + 1 for i, n in enumerate(pn) if n[0] in "NO"]
        peb_idx = [len(pn) + i + 1 for i, n in enumerate(qn) if n[0] in "NO"]
        ndx = f"{OUT}/polar.ndx"
        with open(ndx, "w") as f:
            f.write("[ Prot_polar ]\n" + " ".join(map(str, prot_idx)) + "\n")
            f.write("[ PEB_polar ]\n" + " ".join(map(str, peb_idx)) + "\n")
        run(["mindist", "-s", TPR, "-f", XTC, "-n", ndx, "-d", "0.35", "-dt", "0.1",
             "-od", f"{OUT}/xvg/polar_mindist.xvg", "-on", f"{OUT}/xvg/polar_pairs.xvg"] + tw(),
            stdin_text="0\n1\n", log=f"{OUT}/logs/polar.log")

    # ------------------------------------------------------------ 汇总
    rmsd = read_xvg(f"{OUT}/xvg/rmsd_bb.xvg")
    rg = read_xvg(f"{OUT}/xvg/gyrate.xvg")
    pep = read_xvg(f"{OUT}/xvg/rmsd_peb.xvg")
    mind = read_xvg(f"{OUT}/xvg/mindist_peb.xvg")
    nc = read_xvg(f"{OUT}/xvg/ncontacts.xvg")
    d279 = read_xvg(f"{OUT}/xvg/dist_279_peb.xvg")
    pol = read_xvg(f"{OUT}/xvg/polar_pairs.xvg")

    sc = read_xvg(f"{OUT}/xvg/sasa_complex.xvg") if os.path.exists(f"{OUT}/xvg/sasa_complex.xvg") else None
    sp = read_xvg(f"{OUT}/xvg/sasa_protein.xvg") if os.path.exists(f"{OUT}/xvg/sasa_protein.xvg") else None
    sl = read_xvg(f"{OUT}/xvg/sasa_peb.xvg") if os.path.exists(f"{OUT}/xvg/sasa_peb.xvg") else None
    buried = np.nan
    if sc is not None and sp is not None and sl is not None:
        n = min(len(sc), len(sp), len(sl))
        buried = float(np.mean(sp[:n, 1] + sl[:n, 1] - sc[:n, 1]))

    occ = float((nc[:, 1] > 0).mean())
    pbc_out = int((pep[:, 1] > PBC_CUT).sum())
    pep_med = float(np.median(pep[pep[:, 1] <= PBC_CUT, 1])) if (pep[:, 1] <= PBC_CUT).any() else float("nan")

    row = dict(
        variant=v, atoms=int(n_prot + 93), window_ns=TRIM_TO - TRIM_FROM, ns_done=float(rmsd[-1, 0]),
        rmsd_end_nm=float(rmsd[-1, 1]), rmsd_max_nm=float(rmsd[:, 1].max()),
        Rg_mean_nm=float(rg[:, 1].mean()),
        contact_occupancy=occ, n_contacts_mean=float(nc[:, 1].mean()),
        min_dist_mean_nm=float(mind[:, 1].mean()),
        pep_rmsd_median_nm=pep_med, pbc_outlier_frames=pbc_out,
        d279_PEB_mean_nm=float(d279[:, 1].mean()),
        candidate_hbonds_mean=float(pol[:, 1].mean()),
        buried_sasa_mean_nm2=buried,
    )
    pd.DataFrame([row]).to_csv(f"{OUT}/summary.csv", index=False)

    # 逐变体图
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.4))
    ax[0].plot(rmsd[:, 0], rmsd[:, 1], lw=.7); ax[0].set(xlabel="time (ns)", ylabel="RMSD (nm)", title=f"{v} backbone RMSD")
    ax[1].plot(nc[:, 0], nc[:, 1], lw=.5, color="C3"); ax[1].set(xlabel="time (ns)", ylabel="# contacts < 4 A", title=f"interface contacts (occ {occ*100:.0f}%)")
    ax[2].plot(pol[:, 0], pol[:, 1], lw=.5, color="C6"); ax[2].set(xlabel="time (ns)", ylabel="# polar pairs", title="candidate H-bonds")
    plt.tight_layout(); plt.savefig(f"{OUT}/figs/overview.png", dpi=150); plt.close()
    print(f"      rmsd {row['rmsd_end_nm']:.3f} nm | contacts {row['n_contacts_mean']:.0f} | "
          f"hb {row['candidate_hbonds_mean']:.1f} | buried {buried:.1f} nm2 | pep_med {pep_med:.3f} nm")
    return row


# ---------------------------------------------------------------- 对比
def compare(rows):
    df = pd.DataFrame([r for r in rows if r])
    cmp_dir = f"{AN}/comparison"
    os.makedirs(f"{cmp_dir}/figs", exist_ok=True)
    df.to_csv(f"{cmp_dir}/metrics.csv", index=False)

    ref = df[df.variant == REF]
    if len(ref):
        r0 = ref.iloc[0]
        rank = pd.DataFrame({
            "variant": df.variant,
            "d_contacts_pct": (df.n_contacts_mean / r0.n_contacts_mean - 1) * 100,
            "d_hbonds_pct": (df.candidate_hbonds_mean / r0.candidate_hbonds_mean - 1) * 100,
            "d_buried_pct": (df.buried_sasa_mean_nm2 / r0.buried_sasa_mean_nm2 - 1) * 100,
            "d_pep_shift_nm": df.pep_rmsd_median_nm - r0.pep_rmsd_median_nm,
        })
        # 界面完整性：三项保留率平均（越接近 100% 越"没被破坏"）
        keep = (df.n_contacts_mean / r0.n_contacts_mean +
                df.candidate_hbonds_mean / r0.candidate_hbonds_mean +
                df.buried_sasa_mean_nm2 / r0.buried_sasa_mean_nm2) / 3 * 100
        rank["interface_integrity_pct"] = keep
        rank = rank.sort_values("interface_integrity_pct")
        rank.to_csv(f"{cmp_dir}/ranking.csv", index=False)
        print("\n=== 界面完整性排序（相对 WT，越低=破坏越重）===")
        print(rank.to_string(index=False))

        fig, ax = plt.subplots(1, 3, figsize=(15, 3.8))
        x = np.arange(len(rank))
        ax[0].bar(x, rank.d_contacts_pct, color="C3"); ax[0].set_xticks(x); ax[0].set_xticklabels(rank.variant)
        ax[0].set(ylabel="Δ contacts vs WT (%)", title="界面接触变化")
        ax[1].bar(x, rank.d_hbonds_pct, color="C6"); ax[1].set_xticks(x); ax[1].set_xticklabels(rank.variant)
        ax[1].set(ylabel="Δ H-bonds vs WT (%)", title="候选氢键变化")
        ax[2].bar(x, rank.interface_integrity_pct, color="C0"); ax[2].set_xticks(x); ax[2].set_xticklabels(rank.variant)
        ax[2].axhline(100, ls="--", c="k", lw=.8)
        ax[2].set(ylabel="integrity (%)", title="界面完整性")
        plt.tight_layout(); plt.savefig(f"{cmp_dir}/figs/interface_ranking.png", dpi=150); plt.close()

        # 叠加曲线：接触 / 氢键 / 距离
        fig, ax = plt.subplots(1, 3, figsize=(15, 3.6))
        for _, r in df.iterrows():
            v = r.variant
            nc = read_xvg(f"{AN}/{v}/xvg/ncontacts.xvg")
            ax[0].plot(nc[:, 0], np.convolve(nc[:, 1], np.ones(50) / 50, mode="same"), lw=.9, label=v)
            pol = read_xvg(f"{AN}/{v}/xvg/polar_pairs.xvg")
            ax[1].plot(pol[:, 0], np.convolve(pol[:, 1], np.ones(5) / 5, mode="same"), lw=.9, label=v)
            d = read_xvg(f"{AN}/{v}/xvg/dist_279_peb.xvg")
            ax[2].plot(d[:, 0], d[:, 1], lw=.5, alpha=.7, label=v)
        for a, t, yl in zip(ax, ["interface contacts", "candidate H-bonds", "res279 <-> PEB"],
                            ["# pairs < 4 A", "# polar pairs", "min dist (nm)"]):
            a.set(xlabel="time (ns)", ylabel=yl, title=t); a.legend(fontsize=7)
        plt.tight_layout(); plt.savefig(f"{cmp_dir}/figs/overlay.png", dpi=150); plt.close()
        print(f"\n[cmp] {cmp_dir}/metrics.csv, ranking.csv, figs/")
    return df


def main(argv):
    global TRIM_FROM, TRIM_TO
    reuse = "--reuse" in argv
    args = [a for a in argv if not a.startswith("-")]
    if "--from" in argv:
        TRIM_FROM = float(argv[argv.index("--from") + 1])
    if "--to" in argv:
        TRIM_TO = float(argv[argv.index("--to") + 1])
    # 数值参数不要被当成变体名
    skip = set()
    for flag in ("--from", "--to"):
        if flag in argv:
            skip.add(argv[argv.index(flag) + 1])
    args = [a for a in args if a not in skip]
    variants = args or DEFAULT_VARS
    rows = []
    for v in variants:
        try:
            rows.append(analyze(v, reuse=reuse))
        except Exception as e:
            print(f"[fail] {v}: {str(e)[-400:]}")
            rows.append(None)
    if TRIM_FROM == 0:
        compare(rows)
    else:
        print("[skip] 窗口一致性检验模式，不出跨变体对比（含 --from）")


if __name__ == "__main__":
    main(sys.argv[1:])
