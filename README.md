# StructurePrep-Toolkit

## 📌 项目简介 (Overview)
在计算生物学（如 AutoDock 对接、GROMACS MD 模拟、FoldX 能量计算）中，结构文件的格式缺陷常常导致致命错误。本项目提供了一套轻量级的 Python 处理流程，用于解决真实世界 PDB 数据中的常见坑点，保障下游计算（如稳定性改造、亲和力预测）的准确性。

## ⚠️ 解决的核心痛点 (Problems Solved)
- **配体拓扑丢失（环断裂）**：PyMOL 导出的 `.mol` 文件在进入 AutoDock 或 GROMACS 时经常发生芳香环键级错误，导致对接出错甚至模拟崩溃。本项目使用 ACPYPE 进行了自动化拓扑重构，确保模拟正确。
- **FoldX 忽略配体**：FoldX 严格要求 PDB 格式，包含配体（如 BTN、Kcr）的复合物直接输入常导致配体被静默丢弃。本项目开发了自动化脚本，精准分离蛋白与配体，确保结合能计算的准确性。
- **自动化流程**：将繁琐的命令行手动操作封装为 Python 脚本，减少人为错误，提高蛋白工程改造的效率。

## 📁 项目文件说明 (What's inside)
- `FoldX Ligand Parameterizing.ipynb`：配体参数化实战教程。
- `Foldx_bind_rank.ipynb`：亲和力突变扫描与结构分析实战案例。
- `Rosetta.ipynb`：对突变亲和力进行精细计算。

## 🚀 快速开始 (Quick Start)
*(运行前请确保安装了相关的 Python 环境)*
```bash
pip install -r requirements.txt
```

---

## 📊 结果展示 (Results)

### 🔬 案例 1：CALB 酶与 pNPB 底物的结合模式分析与 MD 灾难排查

#### 1. PyMOL 导出 Mol 格式导致配体芳香环和双键键级丢失等拓扑结构丢失进而导致 MD 模拟灾难
<img width="800" alt="PyMol Results" src="https://github.com/user-attachments/assets/0d93f6f7-d5ce-4296-ba78-d659a243fb1b" />

> **现象与排查**：：在进行 GROMACS MD 模拟前，PyMOL 导出的 `pNPB.mol` 丢失了芳香环和双键的键级（黄色为错误导出，粉色为原始 SDF）。由于拓扑错误导致力场参数异常，配体在 10ns MD 模拟中迅速飘离结合口袋，这正是开发 `StructurePrep-Toolkit` 的核心动机。

#### 2. AutoDock Vina 初始对接结果
<img width="800" alt="Autodock Vina Results" src="https://github.com/user-attachments/assets/b08e5590-a71f-4d34-994f-7881005afa4e" />

> **分析**：AutoDock Vina 预测 pNPB 配体进入 CALB 结合口袋，其**酯键**与催化三联体中的 **SER105** 距离约为 **3.2 Å**，处于理想的催化反应距离。

#### 3. GROMACS 10ns MD 模拟后的真实构象
<img width="800" alt="image" src="https://github.com/user-attachments/assets/0503d39e-3200-4fd9-9dfb-14d9000f2217" />


> **分析**：经过 10ns MD 模拟后，提取距离最小的一帧显示，pNPB 的酯键与 SER105 的距离拉长至 **4.9 Å**，推测 mol 构象错误导致 AutoDock 结果不理想，后续需要继续验证。

> **结论**：PyMOL 导出丢失键级 → 对接输入化学错误 → MD 中底物飘离口袋。


### 🔬 案例 2：hGAS41(20-145) 溶解度工程——突变设计与 MD 验证

> **背景**：由于低溶解度需进行突变设计，疏水斑块分析得到双突变位点 **ILE21Thr & MET128SER**。结合 FoldX 热稳定性分析与 GROMACS 模拟，计算 RMSD 和 RMSF。虽然配体 Kcr 也出现拓扑学问题，但整体依旧结合在口袋，100ns 稳定性较高。

#### 4：FoldX 热稳定性评估
<img width="800" alt="FoldX ddG" src="https://github.com/user-attachments/assets/890b4d75-bf92-4a7f-9450-32b74aae6b4c" />

> **分析**：FoldX 计算 hGAS41(20-145) ILE21Thr & MET128SER 的 \(\Delta\Delta G = 0.221\) kcal/mol，几乎没有破坏蛋白质的整体结构稳定性。

#### 5：结构叠合对比
<img width="800" alt="结构叠合" src="https://github.com/user-attachments/assets/84be34f7-82be-484e-9b1e-3966d776dc66" />

> **说明**：图中 hGAS41(20-145) WT（绿色）与 Kcr（天蓝色）结合晶体解析图，与 hGAS41(20-145) ILE21Thr & MET128SER（粉色）和 Kcr（蓝紫色）GROMACS 模拟图对比。

#### 6：RMSD 分析
<img width="600" alt="RMSD" src="https://github.com/user-attachments/assets/a581fd37-187f-46ba-b35d-93be292a8c2f" />

> **分析**：RMSD 计算 MD 结果稳定，说明双突变可行。

#### 7：RMSF 分析
<img width="600" alt="RMSF" src="https://github.com/user-attachments/assets/3c4b61e4-9263-4773-98a2-c6a156b004de" />

> **分析**：RMSF 柔性区域与 Pymol 对比高度重合。

---

### 🔬 案例 3：hYEATS2-BTN 结合亲和力改造（胞内探针设计）

> **背景**：在 hYEATS2-BTN 结构基础上进行 FoldX PositionScan 饱和突变设计，提升结合亲和力，为胞内探针设计提供参考。

#### FoldX 突变扫描与 Rosetta 交叉验证
#### 8：FoldX 突变扫描结果 (展示 FoldX PositionScan 得到的 ddG 数值，用于评估突变对稳定性和亲和力的影响)
<img width="600" alt="FoldX Results" src="https://github.com/user-attachments/assets/ec3f833c-8080-43fd-bfac-50f1457414c7" />

<br>

#### 9：Rosetta 模拟结果 (展示 Rosetta 计算的能量分数，与 FoldX 结果形成交叉验证)
<img width="500" alt="Rosetta Results" src="https://github.com/user-attachments/assets/deacaf32-366e-4ceb-9729-a3b249355911" />

<br>

#### 10：FoldX/Rosetta 结果汇总 (展示 Rosetta 计算的能量分数，与 FoldX 结果形成交叉验证)
<img width="800" alt="image" src="https://github.com/user-attachments/assets/14cbad28-39be-470e-8451-d13910f616db" />

> **注意**：WT 参照（dG_WT=−13.64 REU），方法误差范围（FoldX ±0.5 / Rosetta ±1 kcal/mol ）
