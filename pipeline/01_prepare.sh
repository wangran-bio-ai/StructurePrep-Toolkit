#!/usr/bin/env bash
# 阶段 1：结构准备 —— 拆链、去水、非标准残基清单、断点/缺原子检查
# 用法: bash 01_prepare.sh input/complex.pdb out/ [CHAIN]
set -euo pipefail
PDB="${1:?need complex.pdb}"
OUT="${2:-out}"
CHAIN="${3:-A}"
GMX="${SP_GMX:-gmx}"
mkdir -p "$OUT"

echo "== 1) 只保留指定链的 ATOM（配体/肽由 itp 提供，不能重复计数）=="
awk -v c="$CHAIN" '/^ATOM/ && substr($0,22,1)==c' "$PDB" > "$OUT/prot.pdb"
echo "   atoms: $(grep -c '^ATOM' "$OUT/prot.pdb")"

echo "== 2) 非标准残基清单 =="
awk '/^(ATOM|HETATM)/ {print substr($0,18,3)}' "$PDB" | sort -u | tr '\n' ' '; echo

echo "== 3) 断点/缺原子（需 gromacs-md 技能的脚本）=="
[ -f "$(dirname "$0")/../scripts/scan_truncated.py" ] && \
  python3 "$(dirname "$0")/../scripts/scan_truncated.py" "$OUT/prot.pdb" || \
  echo "   (scan_truncated.py 未随仓库提供，请从 gromacs-md 技能复制)"

echo "== 4) 进入 pdb2gmx 前，确认无 H、无 altloc、无 HETATM =="
echo "   H 原子: $(grep -c ' H  ' "$OUT/prot.pdb" || true)  HETATM: $(grep -c '^HETATM' "$OUT/prot.pdb" || true)"
