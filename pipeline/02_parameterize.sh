#!/usr/bin/env bash
# 阶段 2：配体参数化 —— ACPYPE(GAFF2 + AM1-BCC)
# 用法: bash 02_parameterize.sh conj.mol2 PEPB 0
#   第 3 个参数是净电荷，必须与 SMILES 形式电荷一致
set -euo pipefail
MOL2="${1:?need mol2}"
NAME="${2:-LIG}"
NETQ="${3:-0}"
WORK="${SP_WORK:-$(pwd)/ligpar}"
ACPY="${SP_ACPY:-acpype}"

mkdir -p "$WORK"; cp -f "$MOL2" "$WORK/"
cd "$WORK"
echo "== acpype: $NAME netq=$NETQ =="
"$ACPY" -i "$(basename "$MOL2")" -c bcc -a gaff2 -n "$NETQ" -b "$NAME"

echo "== 校验 itp =="
ITP="$WORK/${NAME}.acpype/${NAME}_GMX.itp"
[ -f "$ITP" ] || { echo "!! 未生成 $ITP"; exit 1; }
awk '/^\[ *atoms *\]/{f=1;next} /^\[/{f=0} f&&NF>=6{n++;q+=$7} END{printf "   atoms=%d  charge_sum=%+.4f (expect %+d)\n", n, q, '"$NETQ"'}' "$ITP"
