#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
字数统计（Diabetologia 口径）：Sections 1-4。

口径（与 v3.6-v3.12 历代声明值一致）：
  - 取 `## 1. Introduction` 至 `## Declarations`，并在 `## Tables` 截断
  - 剥离：change-log 块引用（行首 `> `）、表格行（`|` 开头）、斜体脚注（`*...*` 整行）、
          标题行（`#`）、`**`/`*` 标记、参考文献编号 `[n]`
单一数据源：manuscript/Manuscript_v*.md

用法：
  python3 scripts/count_words.py manuscript/Manuscript_v3.12_Diabetologia_20261024.md
"""
import re
import sys


def count(path):
    lines = open(path, encoding="utf-8").read().split("\n")

    # 定位正文区间
    start = end = None
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith("## 1. Introduction"):
            start = i
        elif start is not None and s.startswith("## Declarations"):
            end = i
            break
    if start is None or end is None:
        raise SystemExit(f"未找到正文区间：{path}")

    body = lines[start:end]

    kept = []
    for ln in body:
        s = ln.strip()
        if s.startswith(">"):            # change-log 块引用
            continue
        if s.startswith("|"):            # 表格行
            continue
        if s.startswith("#"):            # 标题
            continue
        if s.startswith("*") and s.endswith("*") and len(s) > 2 and "**" not in s[:3]:
            continue                     # 斜体脚注整行
        kept.append(ln)
    text = "\n".join(kept)

    text = re.sub(r"\*\*|\*", "", text)      # 强调标记
    text = re.sub(r"\[[0-9,\s\u2013-]+\]", "", text)  # 引用编号
    words = [w for w in re.split(r"\s+", text) if w.strip()]

    # Research in context 单独计
    ric = None
    for i, ln in enumerate(lines):
        if ln.strip().startswith("## Research in context"):
            j = i + 1
            chunk = []
            while (j < len(lines)
                   and not lines[j].strip().startswith("## ")
                   and not lines[j].strip().startswith("**Keywords**")):
                chunk.append(lines[j])
                j += 1
            ric = len([w for w in re.split(r"\s+", "\n".join(chunk)) if w.strip()])
            break
    return len(words), ric


if __name__ == "__main__":
    for p in sys.argv[1:]:
        n, ric = count(p)
        print(f"{p}\n  Sections 1-4 = {n} words | declared (+70 offset) = {n + 70} | Research in context = {ric}")
