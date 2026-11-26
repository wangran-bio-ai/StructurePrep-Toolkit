# data/

**本目录默认不入库。** 晶体结构属于未发表数据，受隐私/优先权约束，**论文发表前不得上传**。

## 现在能放什么

| 类型 | 现在 | 说明 |
|---|---|---|
| 晶体结构（PDB/mmCIF） | ❌ 禁止 | 未发表，仅在本地/GPU 机器保存 |
| MD 原始数据（xtc/trr/tpr/edr/cpt） | ❌ 禁止 | 体积大，且需要原始结构才能解读 |
| 参数化中间产物（mol2/acpype/sq m） | ❌ 禁止 | 可复现，不必入库 |
| **派生的汇总表**（每变体指标 csv） | ✅ 可以 | 放 `results/`，不放这里 |
| **分析图**（png） | ✅ 可以 | 放 `results/figs/` |

## 论文发表后（解禁流程）

1. 把结构文件放入 `data/structures/`，轨迹放 `data/trajectories/`（建议只放抽帧版）
2. 删除 `.gitignore` 里这两行：
   ```gitignore
   data/**
   !data/README.md
   ```
3. 若只想放部分文件，用白名单逐条放行，例如：
   ```gitignore
   data/**
   !data/README.md
   !data/structures/complex_wt.pdb
   ```
4. 提交前自查：`git status --short | grep -E '\.(pdb|xtc|tpr)$'` 应只剩你确认要公开的那些
5. 推荐给数据配一份 `data/README.md`：来源（PDB ID / 自解）、分辨率、处理步骤

## 数据从哪来

GPU 机器上的原始目录（`GROMACS/<变体>/`）**不要直接 git add**。正确做法：

```bash
# 抽帧 + 居中 + 去周期性跳变，生成适合展示/发布的轻量轨迹
gmx trjconv -s md.tpr -f md_<V>.xtc -o data/trajectories/<V>_10ps.xtc \
            -pbc mol -center -dt 10 -ur compact
# 体积通常从 2.3 GB 降到 150–300 MB
```

## 一句话

**代码与汇总结果随仓库走，原始数据等论文发表后再放。** 这条规则由 `repo-publisher` agent 在每次同步时自动检查。
