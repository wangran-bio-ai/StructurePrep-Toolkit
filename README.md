# StructurePrep-Toolkit

## 📌 项目简介 (Overview)
在计算生物学（如 GROMACS MD 模拟、FoldX 亲和力计算）中，结构文件的格式缺陷常常导致致命错误。本项目提供了一套轻量级的 Python 处理流程，用于解决真实世界 PDB 数据中的常见坑点，保障下游计算（如稳定性改造、亲和力预测）的准确性。

## ⚠️ 解决的核心痛点 (Problems Solved)
1. **配体拓扑丢失（环断裂）**：PyMOL 导出的 `.mol` 文件在进入 GROMACS 时经常发生芳香环键级错误，导致模拟崩溃。本项目使用 RDKit/OpenBabel 进行了自动化拓扑重构。
2. **FoldX 忽略配体**：FoldX 严格要求 PDB 格式，包含配体（如 BTN）的复合物直接输入常导致配体被静默丢弃。本项目开发了自动化脚本，精准分离蛋白与配体，确保结合能计算的准确性。
3. **自动化流程**：将繁琐的命令行手动操作封装为 Python 脚本，减少人为错误，提高蛋白工程改造的效率。

## 📁 项目文件说明 (What's inside)
*   `prep_ligand.py`：修复 mol 文件的拓扑与键级。
*   `prep_complex.py`：分离 PDB 中的蛋白与配体，符合 FoldX 规范。
*   `run_foldx_scan.py`：批量执行 FoldX 突变扫描并提取 CSV 结果。
*   **`analysis_notebook.ipynb`**：👈 **重点！** 包含完整的踩坑复现、修复前后的对比图，以及最终的计算结果分析。

## 🚀 快速开始 (Quick Start)

*(运行前请确保安装了相关的 Python 环境)*
```bash
pip install -r requirements.txt
```

## 📊 结果展示 (Results)

![图 1：PyMOL 导出 Mol 格式导致配体芳香环和双键键级丢失（黄色为 PyMOL 导出，粉色为原始 SDF）。这会导致后续 MD 或 Docking 拓扑构建失败。](https://github.com/user-attachments/assets/0d93f6f7-d5ce-4296-ba78-d659a243fb1b)

<br>

**图 2：FoldX 突变扫描结果** (展示 FoldX PositionScan 得到的 ddG 数值，用于评估突变对稳定性和亲和力的影响)

<img width="685" height="280" alt="FoldX Results" src="https://github.com/user-attachments/assets/ec3f833c-8080-43fd-bfac-50f1457414c7" />

<br>

**图 3：Rosetta 模拟结果** (展示 Rosetta 计算的能量分数，与 FoldX 结果形成交叉验证)

<img width="525" height="195" alt="Rosetta Results" src="https://github.com/user-attachments/assets/deacaf32-366e-4ceb-9729-a3b249355911" />
