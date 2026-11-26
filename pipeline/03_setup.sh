#!/usr/bin/env bash
# 阶段 3：体系搭建 —— pdb2gmx → 拼 itp → 盒子 → 溶剂化 → 离子
# 用法: bash 03_setup.sh out/ prot.pdb [amber99sb-ildn] [1.1] [0.15]
set -euo pipefail
OUT="${1:?need out dir}"
PROT="${2:?need prot.pdb}"
FF="${3:-amber99sb-ildn}"
BOX_D="${4:-1.1}"
CONC="${5:-0.15}"
GMX="${SP_GMX:-gmx}"
mkdir -p "$OUT/mdp" "$OUT/logs_md"; cd "$OUT"; cp -f "../$PROT" prot.pdb 2>/dev/null || cp -f "$PROT" .

$GMX pdb2gmx -f prot.pdb -o prot.gro -p protA.top -ff "$FF" -water tip3p -ignh > logs_md/M3_pdb2gmx.log 2>&1
echo "== pdb2gmx =="; grep -E "Total charge|force field" logs_md/M3_pdb2gmx.log | tail -n 2

echo "== 盒子（dodecahedron 比立方省 ~30% 水）=="
$GMX editconf -f prot.gro -o box.gro -bt dodecahedron -d "$BOX_D" -c > logs_md/M4_editconf.log 2>&1
$GMX solvate -cp box.gro -cs spc216.gro -p protA.top -o solv.gro > logs_md/M5_solvate.log 2>&1
grep -i "solvent molecules" logs_md/M5_solvate.log | tail -n 1

printf 'integrator=steep\nemtol=1000\nemstep=0.01\nnsteps=500\ncutoff-scheme=Verlet\ncoulombtype=PME\nrcoulomb=1.0\nrvdw=1.0\npbc=xyz\n' > mdp/ions.mdp
$GMX grompp -f mdp/ions.mdp -c solv.gro -p protA.top -o ions.tpr -maxwarn 3 > logs_md/M5_grompp_ions.log 2>&1
echo SOL | $GMX genion -s ions.tpr -o solv_ions.gro -p protA.top -pname NA -nname CL -neutral -conc "$CONC" > logs_md/M5_genion.log 2>&1
grep -i "will be added" logs_md/M5_genion.log | head -n 2
echo "== 注意：配体 itp 的 #include 要放在力场 include 之后、第一个 [moleculetype] 之前 =="
