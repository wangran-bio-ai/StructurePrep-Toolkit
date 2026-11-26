# 五阶段流水线入口

每个脚本都是一个薄封装，真正的逻辑在 `notebooks/` 与 `scripts/`。环境变量统一：

```bash
export SP_BASE=/root/autodl-tmp/GROMACS     # 项目根
export SP_VARIANT=R279E                     # 变体名
export SP_GMX=/root/gromacs_cuda/bin/gmx
export SP_NTOMP=8
```

| 脚本 | 阶段 | 做什么 |
|---|---|---|
| `01_prepare.sh` | 结构准备 | 拆链、去水、非标准残基清单、断点/缺原子检查 |
| `02_parameterize.sh` | 配体参数化 | ACPYPE(GAFF2/AM1-BCC) → 校验原子数/电荷 |
| `03_setup.sh` | 体系搭建 | pdb2gmx → 拼 itp → 盒子 → 溶剂化 → 离子 |
| `04_run.sh` | 平衡与生产 | EM → NVT → NPT → 生产（GPU，nohup） |
| `05_analyze.sh` | 界面分析 | analysis_agent 批量跑指标 + 排序 |

用法：

```bash
bash pipeline/01_prepare.sh  input/complex.pdb  out/
bash pipeline/02_parameterize.sh out/conj.mol2 PEPB 0
bash pipeline/03_setup.sh    out/ prot.pdb
bash pipeline/04_run.sh      out/ 50000000     # 步数
bash pipeline/05_analyze.sh  WT Y262A W282A
```

> 每一步的详细说明与判据见 `docs/`，建体系的完整流程见 `notebooks/MD_pipeline.ipynb`。
