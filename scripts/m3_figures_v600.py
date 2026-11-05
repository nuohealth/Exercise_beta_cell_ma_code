#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""m3_figures_v600.py — v6.0 图件：仅重绘 Fig 3（分组条形图）

为什么重绘：v5.0 的 Fig 3 是瀑布图，负的分泌组分用白色填充画在蓝色条之上，
结果视觉上把"分泌下降"渲染成"向右延伸"（负值段的起点在 SI+SEC，终点在 SI，
白底遮住了本应可见的蓝条）。数值坐标是对的，但语义方向反了。
v6.0 改为两条并列的独立条形：SI 与 SEC 各自从 0 起画，负值向左，符号无歧义；
每一臂中点给出观测 Δln DI 的菱形，并在分泌条末端标注符号化数值。

唯一数据源: analysis/arms_v500.json（8 臂；Karstoft 无分泌值，不绘）
输出: figures/Fig3_decomposition_v600.{png,pdf}
用法: python3 scripts/m3_figures_v600.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(BASE, "figures")
VEC = os.path.join(FIG, "vector")
os.makedirs(VEC, exist_ok=True)

CM = 1 / 2.54
PRINT_W = 16.0

plt.rcParams.update({
    "font.family": "Liberation Sans",
    "font.size": 8,
    "axes.linewidth": 0.75,
    "axes.edgecolor": "black",
    "savefig.dpi": 300,
    "figure.dpi": 300,
    "savefig.facecolor": "white",
    "savefig.transparent": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

C_SI = "#1F4E79"
C_SEC = "#E08214"
C_DI = "#000000"

NAME = {"DOSE-EX 2023": "DOSE-EX", "Karstoft 2014": "Karstoft",
        "Slentz 2009": "Slentz", "Lee 2012": "Lee", "Carr 2005": "Carr"}
ARM = {"HED": "HED", "MED": "MED", "LOWMOD": "low/mod", "LOWVIG": "low/vig",
       "HIGHVIG": "high/vig", "IWT": "IWT", "AE": "aerobic", "RE": "resistance", "": ""}
STRAT = {"T2D": "T2D", "NGT": "NGT", "IGT": "IGT"}

ORDER = ["DOSE-EX 2023", "Slentz 2009", "Lee 2012", "Carr 2005"]
ARMORDER = {"DOSE-EX 2023": ["MED", "HED"],
            "Slentz 2009": ["LOWMOD", "LOWVIG", "HIGHVIG"],
            "Lee 2012": ["AE", "RE"], "Carr 2005": [""]}


def save(fig, name):
    p_png = os.path.join(FIG, name + ".png")
    fig.savefig(p_png, dpi=300, facecolor="white", transparent=False)
    fig.savefig(os.path.join(VEC, name + ".pdf"), facecolor="white")
    plt.close(fig)
    im = Image.open(p_png)
    if im.mode != "RGB":
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[-1] if im.mode == "RGBA" else None)
        bg.save(p_png, dpi=(300, 300))
    print(f"  wrote {name}.png + vector pdf")


def load():
    arms = json.load(open(os.path.join(BASE, "analysis", "arms_v500.json"), encoding="utf-8"))
    si = {(r["study"], r["arm"]): r for r in arms["SI"]}
    sec = {(r["study"], r["arm"]): r for r in arms["SEC"]}
    di = {(r["study"], r["arm"]): r for r in arms["DI"]}
    return si, sec, di


def label(study, arm, stratum):
    base = f"{NAME[study]} {ARM[arm]}".strip()
    return f"{base}  [{STRAT[stratum]}]"


def fig3():
    si, sec, di = load()

    labels, sv, cv, dv = [], [], [], []
    for st in ORDER:
        for arm in ARMORDER[st]:
            r = di.get((st, arm))
            if r is None:
                continue
            s_row, c_row = si.get((st, arm)), sec.get((st, arm))
            if s_row is None or c_row is None:
                continue
            labels.append(label(st, arm, r["stratum"]))
            sv.append(s_row["ln"])
            cv.append(c_row["ln"])
            dv.append(r["ln"])

    n = len(labels)
    fig, ax = plt.subplots(figsize=(PRINT_W * CM, (n * 1.00 + 3.1) * CM))
    ys = np.arange(n)[::-1]

    H = 0.34          # bar height
    GAP = 0.19        # half distance between the two bars of one arm

    for y, s, c, d in zip(ys, sv, cv, dv):
        # sensitivity bar (upper), always drawn from zero
        ax.barh(y + GAP, s, height=H, color=C_SI, edgecolor="none", zorder=3)
        # secretion bar (lower), from zero; negative extends left
        ax.barh(y - GAP, c, height=H, color=C_SEC, edgecolor="none", zorder=3)
        # observed total, at the arm's mid-row
        ax.plot(d, y, "D", color=C_DI, ms=4.0, mfc="white", mew=1.05, zorder=6)
        # signed value at the secretion bar's end (makes the sign unmistakable)
        off = 0.012 if c >= 0 else -0.012
        ax.text(c + off, y - GAP, f"{c:+.2f}".replace("-", "\u2212"),
                ha="left" if c >= 0 else "right", va="center",
                fontsize=6.4, color=C_SEC, zorder=7)
        # sensitivity value, kept lighter to avoid clutter
        ax.text(s + 0.012, y + GAP, f"{s:+.2f}".replace("-", "\u2212"),
                ha="left", va="center", fontsize=6.4, color=C_SI, zorder=7)

    ax.axvline(0, color="#404040", lw=0.8, zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=7.2)
    ax.set_ylim(-0.85, n - 0.15)
    ax.set_xlim(-0.30, 0.98)
    ax.set_xlabel("Change vs control (\u0394ln unit)", fontsize=7.8, labelpad=3)
    ax.tick_params(axis="x", width=0.75, length=2.5, labelsize=7.4)
    ax.tick_params(axis="y", length=0)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

    handles = [Patch(facecolor=C_SI, label="insulin sensitivity"),
               Patch(facecolor=C_SEC, label="insulin secretion"),
               Line2D([0], [0], marker="D", color="none", mfc="white", mec="black",
                      mew=1.05, ms=4.0, label="observed \u0394ln DI")]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.45, 1.005),
              ncol=3, frameon=False, fontsize=7.0, handlelength=1.3,
              columnspacing=1.5, handletextpad=0.5)

    fig.subplots_adjust(left=0.275, right=0.985, top=0.878, bottom=0.115)
    save(fig, "Fig3_decomposition_v600")


if __name__ == "__main__":
    fig3()
