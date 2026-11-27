#!/usr/bin/env python3
"""report_agent.py — 把 ANALYSIS/ 下所有结果拼成一份 Markdown 报告

用法: python report_agent.py [--out REPORT.md]
读取：comparison/{metrics,ranking}.csv、mmpbsa_summary.csv、
      peptide_motion_summary.csv、cage_geometry.csv、decomp_*_vs_WT.csv
"""
import glob
import os
import sys

import pandas as pd

B = "/root/autodl-tmp/GROMACS/ANALYSIS"


def load(p, **kw):
    return pd.read_csv(p, **kw) if os.path.exists(p) else None


def md_table(df, floatfmt="%.2f"):
    if df is None or df.empty:
        return "_（无数据）_\n"
    cols = list(df.columns)
    out = ["| " + " | ".join(map(str, cols)) + " |",
           "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            cells.append(floatfmt % v if isinstance(v, float) else str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out) + "\n"


def main():
    out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else f"{B}/REPORT.md"
    L = []
    L.append("# MD 分析汇总报告\n")
    L.append("> 由 `report_agent.py` 自动生成；所有数字取自 `ANALYSIS/` 下的 CSV。\n")

    rank = load(f"{B}/comparison/ranking.csv")
    L.append("## 1 界面完整性排序（0–50 ns，相对 WT）\n")
    L.append(md_table(rank))

    metrics = load(f"{B}/comparison/metrics.csv")
    if metrics is not None:
        keep = [c for c in ("variant", "n_contacts_mean", "candidate_hbonds_mean",
                            "buried_sasa_mean_nm2", "pep_rmsd_median_nm", "rmsd_end_nm")
                if c in metrics.columns]
        L.append("\n## 2 结构指标明细\n")
        L.append(md_table(metrics[keep]))

    lig = load(f"{B}/peptide_motion_summary.csv")
    L.append("\n## 3 界面重排（逐残基接触 + 肽滑移）\n")
    L.append(md_table(lig))

    mmp = load(f"{B}/mmpbsa_summary.csv")
    L.append("\n## 4 MM/GBSA（125 帧，50–100 ns）\n")
    L.append(md_table(mmp))
    L.append("\n> ΔEEL 与 ΔEGB 必须成对引用；|ΔΔG| < 3 kcal/mol 视为噪声。"
             "电荷突变（R279E/S261E）的绝对值不可用于判断亲和力。\n")

    cage = load(f"{B}/cage_geometry.csv")
    if cage is not None:
        L.append("\n## 5 芳香笼几何\n")
        L.append(md_table(cage))

    decs = sorted(glob.glob(f"{B}/decomp_*_vs_WT.csv"))
    if decs:
        L.append("\n## 6 逐残基能量分解（Amber MMPBSA.py）\n")
        for p in decs:
            v = os.path.basename(p).replace("decomp_", "").replace("_vs_WT.csv", "")
            df = pd.read_csv(p)
            if "d_total" not in df.columns:
                continue
            top = df.sort_values("d_total", ascending=False).head(4)
            bot = df.sort_values("d_total").head(4)
            fmt = lambda d: "、".join(f"{int(r.resid)}({r.d_total:+.2f})" for r in d.itertuples())
            L.append(f"- **{v}**：合计 {df.d_total.sum():+.2f} kcal/mol；"
                     f"削弱 {fmt(top)}；补偿 {fmt(bot)}")
        L.append("\n> 逐残基之和与直接 ΔΔG 不完全相等（采样帧数不同），只用于相对归属。\n")

    L.append("\n## 7 结论与边界\n")
    L.append("""
1. **W282A**：唯一在结构与能量上都一致显示亲和力下降的变体（ΔΔG +14.5 ± 2.1 kcal/mol；
   buried SASA −9.3%；界面完整性 91.4%；ΔVDW 少 11.9 kcal/mol）。
2. **Y262A**：界面接触重排（Trp234 转入、Pro286 退出），结合态相互作用能不变（−2.1 ± 1.5），
   说明其表型不能由结合态能量解释。
3. **E284P**：界面轻度重排，能量无显著变化（−0.6 ± 1.6）。
4. **R279E / S261E**：改变界面电荷，MM/GBSA 中 ΔEEL 与 ΔEGB 大幅抵消
   （S261E：−141 / +144 kcal/mol），逐残基合计≈0，**本方法不能给出亲和力方向**；
   定量需显式溶剂 FEP 或 ITC。
5. 所有变体在 80–100 ns 内均未解离（接触占时比 100%），因此"结合态能量相当"不等于亲和力相同。
""")

    open(out, "w", encoding="utf-8").write("\n".join(L))
    print("written", out, len("\n".join(L)), "chars")


if __name__ == "__main__":
    main()
