# results/

**只放汇总表与图，不放原始数据。**

| 文件 | 内容 |
|---|---|
| `metrics.csv` | 每个变体的全指标（含 `window_ns`） |
| `ranking.csv` | 界面完整性排序（相对参考变体归一化） |
| `figs/interface_ranking.png` | 接触变化 / 氢键变化 / 完整性 三柱图 |
| `figs/overlay.png` | 多变体接触、氢键、关键距离叠加曲线 |

## 当前状态

三变体（WT / Y262A / W282A），统一 0–50 ns 窗口；Y262A 与 W282A 的 100 ns 轨迹与
E284P、S261E 的数据在补，更新后会替换本目录下的 csv 与图。

## 更新方式

```bash
python ../scripts/analysis_agent.py WT Y262A W282A E284P S261E --ref WT
cp ANALYSIS/comparison/{metrics,ranking}.csv results/
cp ANALYSIS/comparison/figs/*.png results/figs/
```

> 只提交 `metrics.csv` / `ranking.csv` / `*.png`；`.gitignore` 会拦掉轨迹与结构文件。
