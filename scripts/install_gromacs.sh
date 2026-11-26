#!/usr/bin/env bash
# install_gromacs.sh — 把 GROMACS 装到数据盘（默认 CPU 版，conda-forge/清华源）
#
#   bash install_gromacs.sh                      # → /root/autodl-tmp/md
#   bash install_gromacs.sh /root/autodl-tmp/md 2024.4   # 指定版本
#
# 只在用户明确说"安装"时运行。副作用：新建 $PREFIX（约 1–2 GB），不碰任何已有目录。
# 卸载：rm -rf $PREFIX   （只毁掉这个环境的 gmx 和包缓存）
#
# 装完必须显式传线程数：nproc 在容器里常返回宿主核数（96），实际只给 12 核。
#   gmx mdrun ... -ntmpi 1 -ntomp 12
set -euo pipefail

PREFIX="${1:-/root/autodl-tmp/md}"
VER="${2:-}"
ROOT="$(dirname "$PREFIX")"

echo "=== 侦察 ==="
uname -sr || true
echo "nproc(可能是宿主核数): $(nproc)"
df -h "$ROOT" | tail -1
command -v gmx || echo "gmx: MISSING"
command -v nvcc || echo "nvcc: MISSING（→ 只能先装 CPU 版）"
nvidia-smi --query-gpu=name,driver_version,compute_cap --format=csv,noheader 2>/dev/null || echo "无 GPU"

if [ -x "$PREFIX/bin/gmx" ]; then
  echo "=== 已存在，不覆盖：$PREFIX/bin/gmx 已可执行（要重装先 rm -rf $PREFIX）==="
  "$PREFIX/bin/gmx" --version | head -3
  exit 0
fi

echo "=== micromamba ==="
mkdir -p "$ROOT/bin"
if [ ! -x "$ROOT/bin/micromamba" ]; then
  curl -Ls --max-time 180 https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$ROOT" bin/micromamba
fi

PKG="gromacs"
[ -n "$VER" ] && PKG="gromacs=$VER"

echo "=== 安装 $PKG → $PREFIX ==="
"$ROOT/bin/micromamba" create -y -p "$PREFIX" \
  -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge \
  -c conda-forge ${PKG} python=3.11 2>&1 | tail -15

echo "=== 结果 ==="
"$PREFIX/bin/gmx" --version 2>&1 | head -8
df -h "$ROOT" | tail -1
echo "提示：conda-forge 版是 CPU-only（无 -nb gpu），生产前需源码编译 CUDA 版。"
