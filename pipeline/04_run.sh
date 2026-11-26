#!/usr/bin/env bash
# 阶段 4：EM → NVT → NPT → 生产
# 用法: bash 04_run.sh out/ [生产步数=50000000] [ntomp=8]
set -euo pipefail
OUT="${1:?need out dir}"
NSTEPS="${2:-50000000}"
NTOMP="${3:-8}"
GMX="${SP_GMX:-gmx}"
cd "$OUT"

cat > mdp/em.mdp <<'EOF'
integrator=steep
emtol=1000
emstep=0.01
nsteps=50000
cutoff-scheme=Verlet
coulombtype=PME
rcoulomb=1.0
rvdw=1.0
pbc=xyz
EOF
cp mdp/em.mdp mdp/ions.mdp

$GMX grompp -f mdp/em.mdp -c solv_ions.gro -p protA.top -o em.tpr -maxwarn 5 > logs_md/M5c_grompp_em.log 2>&1
$GMX mdrun -deffnm em -ntmpi 1 -ntomp "$NTOMP" -nb gpu -pin off > logs_md/M5d_mdrun_em.log 2>&1
grep -E "converged|Maximum force" logs_md/M5d_mdrun_em.log | tail -n 2

cat > mdp/nvt.mdp <<'EOF'
integrator=md
dt=0.002
nsteps=50000
constraint_algorithm=lincs
constraints=h-bonds
lincs_iter=1
lincs_order=4
cutoff-scheme=Verlet
coulombtype=PME
rcoulomb=1.0
rvdw=1.0
DispCorr=EnerPres
tcoupl=V-rescale
tc-grps=Protein Non-Protein
ref_t=300 300
tau_t=0.1 0.1
pcoupl=no
pbc=xyz
gen_vel=yes
gen_temp=300
gen_seed=-1
EOF
sed -e 's/^gen_vel=yes/gen_vel=no/' -e 's/^pcoupl=no/pcoupl=C-rescale\npcoupltype=isotropic\ntau_p=2.0\nref_p=1.0\ncompressibility=4.5e-5/' mdp/nvt.mdp > mdp/npt.mdp
sed -e "s/^nsteps=50000/nsteps=$NSTEPS/" mdp/npt.mdp > mdp/md.mdp
printf 'nstxout-compressed=5000\ncompressed-x-grps=System\nnstenergy=5000\nnstlog=5000\n' >> mdp/md.mdp

for pair in "nvt em.gro none" "npt nvt.gro nvt.cpt"; do
  set -- $pair
  if [ "$3" = "none" ]; then
    $GMX grompp -f mdp/$1.mdp -c "$2" -p protA.top -o $1.tpr -maxwarn 5 > logs_md/M_$1.log 2>&1
  else
    $GMX grompp -f mdp/$1.mdp -c "$2" -t "$3" -p protA.top -o $1.tpr -maxwarn 5 > logs_md/M_$1.log 2>&1
  fi
  $GMX mdrun -deffnm $1 -ntmpi 1 -ntomp "$NTOMP" -nb gpu -pme gpu -pin off > logs_md/M${1}_run.log 2>&1
  grep -m1 Performance logs_md/M${1}_run.log || true
done

$GMX grompp -f mdp/md.mdp -c npt.gro -t npt.cpt -p protA.top -o md.tpr -maxwarn 5 > logs_md/M7_grompp.log 2>&1
echo "== 生产（后台）=="
nohup $GMX mdrun -deffnm md -s md.tpr -ntmpi 1 -ntomp "$NTOMP" -nb gpu -pme gpu -pin off > logs_md/M7_mdrun.log 2>&1 &
echo "  pid=$!  tpr=md.tpr  输出=md.*（-deffnm 只改输出名）"
