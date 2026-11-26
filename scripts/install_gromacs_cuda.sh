#!/usr/bin/env bash
# install_gromacs_cuda.sh — 源码编译带 CUDA 的 GROMACS（生产用）
#
#   bash install_gromacs_cuda.sh            # 默认 2024.4 -> /root/gromacs_cuda_244
#   bash install_gromacs_cuda.sh 2024.4 /root/gromacs_cuda_244
#
# 关键：-DGMX_SIMD=AVX2_256 必须写死。用 -march=native 编出来的 AVX-512 二进制
# 换到只有 AVX2 的机器上会 Illegal instruction (core dumped) / rc=132。
# 只在用户明确说"编译/安装"时运行；约 8–25 min（12 核）。
set -euo pipefail

VER="${1:-2024.4}"
PREFIX="${2:-/root/gromacs_cuda_$(echo $VER | tr -d '.')}"
SRCROOT=/root/autodl-tmp/src
SRC=$SRCROOT/gromacs-$VER
BUILD=$SRC/build
SM="${GMX_CUDA_SM:-86}"     # 3090 = 86

echo "=== preflight ==="
gcc --version 2>/dev/null | head -n 1 || echo "gcc MISSING"
cmake --version 2>/dev/null | head -n 1 || echo "cmake MISSING"
/usr/local/cuda/bin/nvcc --version 2>/dev/null | tail -n 2 || echo "nvcc MISSING"
df -h "$SRCROOT" | tail -n 1

if ! command -v gcc >/dev/null 2>&1; then
  apt-get update -qq || true; apt-get install -y -qq build-essential
fi
if ! command -v cmake >/dev/null 2>&1; then
  /root/autodl-tmp/bin/micromamba create -y -p /root/autodl-tmp/build -c conda-forge cmake make
  export PATH=/root/autodl-tmp/build/bin:$PATH
fi
# FFTW：优先系统包（快），否则让 GROMACS 自建（下载慢）
OWN_FFTW=ON
if ! ls /usr/include/fftw3.h >/dev/null 2>&1; then
  apt-get update -qq || true
  apt-get install -y -qq libfftw3-dev 2>/dev/null || true
fi
[ -f /usr/include/fftw3.h ] && OWN_FFTW=OFF

echo "=== source ==="
mkdir -p "$SRCROOT"
[ -f "$SRC/CMakeLists.txt" ] || {
  cd "$SRCROOT"
  for U in https://ftp.gromacs.org/gromacs/gromacs-$VER.tar.gz \
           https://mirrors.tuna.tsinghua.edu.cn/gromacs/gromacs-$VER.tar.gz \
           https://github.com/gromacs/gromacs/archive/refs/tags/v$VER.tar.gz ; do
    echo "-- $U"
    curl -fL --max-time 300 -o "$SRCROOT/gmx.tar.gz" "$U" && break || true
  done
  tar xzf "$SRCROOT/gmx.tar.gz" -C "$SRCROOT"
  [ -f "$SRC/CMakeLists.txt" ] || { echo "!! source missing"; exit 1; }
}

echo "=== configure ==="
mkdir -p "$BUILD"; cd "$BUILD"
cmake "$SRC" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DGMX_GPU=CUDA \
  -DGMX_SIMD=AVX2_256 \
  -DGMX_CUDA_TARGET_SM="$SM" \
  -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc \
  -DGMX_BUILD_OWN_FFTW="$OWN_FFTW" \
  -DREGRESSIONTEST_DOWNLOAD=OFF \
  -DGMX_BUILD_UNITTESTS=OFF \
  -DBUILD_TESTING=OFF \
  -DGMX_DEFAULT_SUFFIX=OFF \
  -DGMX_X11=OFF \
  > /root/cmake_$VER.log 2>&1 || { echo "!! cmake failed"; tail -n 30 /root/cmake_$VER.log; exit 1; }
grep -iE "SIMD instructions|GPU support" /root/cmake_$VER.log | head -n 3 || true

echo "=== make -j$(nproc | xargs -I{} sh -c 'echo 12') ==="
make -j12 > /root/make_$VER.log 2>&1 || { echo "!! make failed"; tail -n 30 /root/make_$VER.log; exit 1; }
make install > /root/install_$VER.log 2>&1 || { echo "!! install failed"; tail -n 20 /root/install_$VER.log; exit 1; }

echo "=== verify ==="
"$PREFIX/bin/gmx" --version 2>&1 | head -n 8
"$PREFIX/bin/gmx" --version 2>&1 | grep -iE "gpu|simd" | head -n 4
echo "add to PATH:  echo 'export PATH=$PREFIX/bin:\$PATH' >> ~/.bashrc"
