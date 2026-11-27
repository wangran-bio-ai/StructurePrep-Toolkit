#!/usr/bin/env python3
"""spconfig.py — 统一配置加载（阈值/权重/窗口），避免经验数值硬编码在公开代码里。

查找顺序（先到先用）：
  1. 环境变量 SP_CONFIG 指定的文件
  2. ./analysis.ini          （项目/运行目录，已被 .gitignore 忽略）
  3. ../analysis.ini         （仓库上一级）
  4. ~/.structureprep/analysis.ini
  5. 仓库自带 analysis.example.ini（**通用占位值，不是校准值**）

用法：
    from spconfig import cfg
    cut = cfg.getfloat("contacts", "cutoff_nm", 0.4)
    w   = cfg.getfloat("integrity", "weight_contacts", 1.0)

公开仓库里只保留 analysis.example.ini 的结构；校准值放在私有配置中。
"""
import configparser
import os

_BASES = [
    os.environ.get("SP_CONFIG", ""),
    "analysis.ini",
    os.path.join("..", "analysis.ini"),
    os.path.expanduser("~/.structureprep/analysis.ini"),
]
_FALLBACK = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis.example.ini")


class Config:
    def __init__(self, path=None):
        self.cp = configparser.ConfigParser()
        self.source = None
        for cand in ([path] if path else []) + _BASES:
            if cand and os.path.exists(cand):
                self.cp.read(cand, encoding="utf-8")
                self.source = os.path.abspath(cand)
                break
        if self.source is None and os.path.exists(_FALLBACK):
            self.cp.read(_FALLBACK, encoding="utf-8")
            self.source = f"{os.path.abspath(_FALLBACK)} (example defaults)"

    def get(self, sec, key, default=None):
        return self.cp.get(sec, key, fallback=default)

    def getfloat(self, sec, key, default=0.0):
        return self.cp.getfloat(sec, key, fallback=default)

    def getint(self, sec, key, default=0):
        return self.cp.getint(sec, key, fallback=default)

    def sections(self):
        return self.cp.sections()


cfg = Config()

if __name__ == "__main__":
    print("config source:", cfg.source)
    for s in cfg.sections():
        for k, v in cfg.cp[s].items():
            print(f"  [{s}] {k} = {v}")
