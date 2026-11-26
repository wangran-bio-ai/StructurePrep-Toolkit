# StructurePrep-Toolkit

面向**蛋白–多肽/配体复合物**的结构准备与分子动力学分析流水线。GROMACS 2024.4 + Amber 力场族（amber99sb-ildn / tip3p / GAFF2）。

设计目标：把"拿到一个晶体结构 → 能跑的 MD 体系 → 可解释的界面指标"这条链路做成**可复现、可审计**的五个阶段，每一步都有明确的通过判据。

---

## Pipeline

```
01 结构准备      晶体 PDB → 拆链/去水/补全 → 可进 pdb2gmx 的坐标
02 配体参数化    配体/共价修饰肽 → GAFF2 + AM1-BCC → 可直接并入的 itp
03 体系搭建      pdb2gmx → 拼 itp → 盒子 → 溶剂化 → 离子
04 平衡与生产    EM → NVT → NPT → 生产（GPU）
05 界面分析      RMSD/Rg/RMSF、接触、buried SASA、候选氢键、完整性排序
```

| 阶段 | 入口 | 通过判据 |
|---|---|---|
| 01 结构准备 | `pipeline/01_prepare.sh` | 0 断点、0 缺原子、非标准残基清单明确 |
| 02 配体参数化 | `notebooks/ligand_param.ipynb` | 原子数、电荷和 = 净电荷、residue 名一致 |
| 03 体系搭建 | `notebooks/MD_pipeline.ipynb`（M3） | `Total charge` 对账、水数 ×3 + 溶质 = gro 原子数 |
| 04 平衡与生产 | 同上（M4–M7） | EM `Fmax < 1000`、NVT ±5 K、NPT 密度 ≈1000 kg/m³ |
| 05 界面分析 | `scripts/analysis_agent.py` | 接触占时比、buried SASA、候选氢键、完整性排序 |

---

## 快速开始

```bash
# 依赖
micromamba create -p ./env -c conda-forge gromacs=2024.4 ambertools acpype openbabel python=3.11
# 或源码编译 CUDA 版
bash scripts/install_gromacs_cuda.sh

# 05 分析：一条命令跑完所有变体
python scripts/analysis_agent.py WT Y262A W282A --ref WT
```

`analysis_agent.py` 的输出：

```
ANALYSIS/
├── <V>/{summary.csv, xvg/, figs/, logs/}
└── comparison/{metrics.csv, ranking.csv, figs/{interface_ranking,overlay}.png}
```

`ranking.csv` 里的 `interface_integrity_pct` **越低 = 界面破坏越重**。

---

## 已有结果（蛋白–生物素化肽，五个变体）

见 `results/`。当前口径：**同一盒子、同一力场、统一 0–50 ns 分析窗口**。

<!-- RESULTS:START -->
| variant | window_ns | n_contacts_mean | candidate_hbonds_mean | buried_sasa_mean_nm2 | pep_rmsd_median_nm |
|---|---|---|---|---|---|
| WT | 50 | 480.0 | 9.28 | 10.62 | 0.463 |
| Y262A | 50 | 396.9 | 9.11 | 10.18 | 0.377 |
| W282A | 50 | 421.0 | 8.89 | 9.64 | 0.293 |
<!-- RESULTS:END -->

> ⚠ 结构指标只能说明"界面被削弱"，**不能说明亲和力强弱**。50 ns 普通 MD 看不到 ms 级的解离事件；结合强弱需 ITC 或 FEP/PMF 判定。详见 `docs/metrics_guide.md`。

---

## 数据政策

本仓库**只放代码、文档与汇总结果**：

- ❌ 不上传晶体结构（PDB/mmCIF）、轨迹（xtc/trr）、拓扑（tpr/gro/top）、能量文件（edr/cpt）
- ✅ 只上传：脚本、notebook、汇总表（csv）、图（png）、配置样例

原因：晶体坐标可能未公开，轨迹体积过大。`.gitignore` 已强制拦截。

---

## 目录

```
├── data/            ★ 默认不入库：晶体结构与轨迹（论文发表后再放，见 data/README.md）
├── docs/           方法学、指标判读、常见坑
├── pipeline/       五阶段入口脚本
├── notebooks/      MD_pipeline / ligand_param / Analysis_notebook
├── scripts/        analysis_agent、hbond_geometric、环境安装
└── results/        汇总表与图（无原始数据）
```

> **`data/` 现在是空的占位目录**：晶体结构属未发表数据，论文发表前不公开。解禁步骤见 `data/README.md`。

## 引用与致谢

GROMACS、AmberTools、ACPYPE、OpenBabel。若使用本流水线，请注明来源。
