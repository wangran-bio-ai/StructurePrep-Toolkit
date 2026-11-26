#!/usr/bin/env bash
# 阶段 5：界面分析
# 用法: bash 05_analyze.sh WT Y262A W282A [--ref WT]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
python3 "$HERE/../scripts/analysis_agent.py" "$@"
