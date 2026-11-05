#!/usr/bin/env python3
"""
Fig2：证据缺口图（Evidence Gap Map）
--------------------------------------
用途：把"我们能证明什么、不能证明什么"画成一张不可争辩的矩阵。
每一格 = 该设计 × 该测量组合下已发表研究的数量（N=0 的格子就是缺口）。
数据来源：2026-10-21 PubMed/EuropePMC 检索，PMID 已逐个核对标题与摘要。
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle

for cand in ["Noto Sans CJK SC", "Noto Sans CJK JP", "WenQuanYi Zen Hei", "WenQuanYi Micro Hei"]:
    if any(cand in f.name for f in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = cand
        break
plt.rcParams["axes.unicode_minus"] = False

rows = [
    "L1 运动 → 分泌/敏感性分解\n（人体试验）",
    "L2 Metrnl → β 细胞\n（机制）",
    "L3 运动 → 循环 Metrnl\n（人体干预）",
    "L4 运动 → 循环 Metrnl → β 细胞功能\n（三层串联）",
]
cols = ["横断面\n相关", "单组前后\n（无对照）", "RCT\n（组均值）", "RCT + 动态指标\n（钳夹/IVGTT）", "RCT + 动态 + 可分解\n（几何均值）"]

# 计数
N = np.array([
    [12, 8, 26, 10, 2],   # L1
    [1, 0, 0, 5, 0],      # L2
    [0, 3, 10, 0, 0],     # L3
    [0, 0, 0, 0, 0],      # L4
], dtype=float)

note = [
    ["流行病学", "单臂训练", "26 项\n双报告", "钳夹/IVGTT\n10 项", "Slentz09\nDOSE-EX\n仅 2 项"],
    ["阴性\n(HOMA-β)", "", "", "体外/动物\n5 篇", "—"],
    ["", "无对照\n3 篇", "10 项\n均未测 DI", "", "全空"],
    ["", "", "", "", "全球 0 篇"],
]

fig, ax = plt.subplots(figsize=(12.6, 6.4))
cmap = plt.get_cmap("YlOrRd")
vmax = 26
for i in range(N.shape[0]):
    for j in range(N.shape[1]):
        v = N[i, j]
        if v == 0:
            fc, ec, lw = "#F2F4F6", "#9AA5B1", 1.5
            if note[i][j]:
                fc, ec, lw = "#FDECEA", "#C0392B", 2.2
        else:
            fc = cmap(min(v / vmax, 1.0) * 0.72 + 0.12)
            ec, lw = "white", 1.0
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=fc, edgecolor=ec, linewidth=lw))
        dark = v >= 8
        txt = note[i][j] if note[i][j] else ("0" if v == 0 else "")
        ax.text(j + 0.5, i + 0.44 if v > 0 else i + 0.5, txt, ha="center", va="center",
                fontsize=7.6 if len(txt) > 5 else 10.5,
                color=("white" if dark else "#7B241C") if v == 0 or dark else ("#5D6D7E" if v == 0 else "#1B2631"),
                fontweight="bold" if v > 0 and len(txt) <= 4 else "normal")
        if v > 0:
            ax.text(j + 0.5, i + 0.80, f"n = {int(v)}", ha="center", va="center",
                    fontsize=8.0, color="white" if dark else "#1B2631", fontweight="bold")

ax.set_xlim(0, len(cols)); ax.set_ylim(0, len(rows))
ax.set_xticks(np.arange(len(cols)) + 0.5)
ax.set_xticklabels(cols, fontsize=8.6)
ax.set_yticks(np.arange(len(rows)) + 0.5)
ax.set_yticklabels(rows, fontsize=8.8)
ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)

# 标出三个红格（本文要填的洞）
ax.add_patch(Rectangle((0, 3), 5, 1, facecolor="none", edgecolor="#C0392B", linewidth=3.2, zorder=6))
ax.text(5.08, 3.5, "本文\n主攻", ha="left", va="center", fontsize=9.2, color="#C0392B", fontweight="bold")
ax.add_patch(Rectangle((3, 2), 2, 1, facecolor="none", edgecolor="#E67E22", linewidth=2.4, ls="--", zorder=6))
ax.text(5.08, 2.5, "本文\n顺带填", ha="left", va="center", fontsize=8.6, color="#CA6F1E")
ax.add_patch(Rectangle((4, 0), 1, 1, facecolor="none", edgecolor="#1F6FB4", linewidth=2.4, ls="--", zorder=6))
ax.text(5.08, 0.5, "本文\n核心方法", ha="left", va="center", fontsize=8.6, color="#1F6FB4")

ax.set_title("证据缺口图：四层证据，只有底层是满的，顶层完全空白",
             fontsize=13, loc="left", fontweight="bold", pad=34)
ax.text(0, -0.62, "红框 = 完全空白（本文要填）；橙框 = 部分空白；n = 该格已发表研究数\n"
                  "同行比对可见：运动→Metrnl 有人测，Metrnl→β 细胞有人测，但『运动→Metrnl→DI 分解』三层串联全球 0 篇",
        fontsize=8.2, color="#4D5656", va="top")

plt.tight_layout()
out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")
os.makedirs(out, exist_ok=True)
p = os.path.join(out, "Fig2_evidence_gap_map.png")
plt.savefig(p, dpi=320, bbox_inches="tight", facecolor="white")
print("[figure]", p)
for i, r in enumerate(rows):
    print(f"{r.replace(chr(10),' '):46s}", " ".join(f"{int(x):3d}" for x in N[i]))
