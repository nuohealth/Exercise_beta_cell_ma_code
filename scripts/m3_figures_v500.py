#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""m3_figures_v500.py — v5.0 图件（9 臂 / 5 研究 / 三层）

唯一数据源: analysis/arms_v500.json（由 scripts/arms_v500.py 生成并自校验）
合并估计:   analysis/meta_A_a1p.json

与 v313 的差别只有三处（画布规范、配色、明度分离全部沿用）：
  1. 新增第五个研究 Carr 2005（IGT 层，单臂，无研究内平行臂）；
  2. 分层由两层（T2D / NGT）变三层（T2D / NGT / IGT）；
  3. k 由 8 臂变 9 臂（分泌轴仍 8 臂 —— Karstoft 无分泌值，Carr 补入）。

用法: python3 scripts/m3_figures_v500.py
输出: figures/Fig2_forest_three_axes_v500.{png,pdf}
      figures/Fig3_decomposition_v500.{png,pdf}
"""
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory
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

ORDER = ["DOSE-EX 2023", "Karstoft 2014", "Slentz 2009", "Lee 2012", "Carr 2005"]
ARMORDER = {"DOSE-EX 2023": ["MED", "HED"], "Karstoft 2014": ["IWT"],
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
    pooled = json.load(open(os.path.join(BASE, "analysis", "meta_A_a1p.json"), encoding="utf-8"))
    return arms, pooled


def label(study, arm, stratum):
    base = f"{NAME[study]} {ARM[arm]}".strip()
    return f"{base}  [{STRAT[stratum]}]"


# ══════════════════════════════════════════════════════════════════
# FIG 2 — 三轴森林图（DI 9 臂 / SI 9 臂 / SEC 8 臂）
# ══════════════════════════════════════════════════════════════════
def fig2():
    arms, pooled = load()
    sections = [("DI", "Disposition index", C_DI, "a"),
                ("SI", "Insulin sensitivity", C_SI, "b"),
                ("SEC", "Insulin secretion", C_SEC, "c")]

    ROW, HDR_GAP, POOL_GAP = 0.33, 0.85, 1.05
    rows, y = [], 0.0

    for q, title, col, letter in sections:
        rows.append((y, "header", (title, col, letter)))
        y -= HDR_GAP
        for st in ORDER:
            for arm in ARMORDER[st]:
                r = next((x for x in arms[q] if x["study"] == st and x["arm"] == arm), None)
                if r is None:
                    continue
                lo = (math.exp(r["ln"] - 1.959964 * r["se"]) - 1) * 100
                hi = (math.exp(r["ln"] + 1.959964 * r["se"]) - 1) * 100
                rows.append((y, "arm", (label(st, arm, r["stratum"]),
                                        (math.exp(r["ln"]) - 1) * 100, lo, hi, col)))
                y -= ROW
        m = pooled[q]["new"]
        rows.append((y, "pooled", (f"Pooled (k={m['k']})", m["pct"], m["ci"][0], m["ci"][1], col)))
        y -= POOL_GAP

    DEPTH = -y
    fig, ax = plt.subplots(figsize=(PRINT_W * CM, (DEPTH + 1.1) * CM))
    tr = blended_transform_factory(ax.transAxes, ax.transData)

    for (yy, kind, payload) in rows:
        if kind == "header":
            title, col, letter = payload
            ax.text(-0.012, yy, title, ha="right", va="center", transform=tr,
                    fontsize=8.4, fontweight="bold", color=col, clip_on=False)
            ax.text(-0.012, yy + 0.34, letter, ha="right", va="center", transform=tr,
                    fontsize=8.4, fontweight="bold", color="black", clip_on=False)
            ax.axhline(yy - 0.50, color=col, lw=0.75, alpha=0.45, xmin=0, xmax=1)
        elif kind == "arm":
            lab, p, lo, hi, col = payload
            ax.plot([lo, hi], [yy, yy], color=col, lw=1.0, solid_capstyle="butt", zorder=2)
            ax.plot(p, yy, "s", color=col, ms=3.6, zorder=3)
            ax.text(-0.012, yy, lab, ha="right", va="center", transform=tr,
                    fontsize=7.2, clip_on=False)
        else:
            lab, p, lo, hi, col = payload
            ax.plot([lo, hi], [yy, yy], color=col, lw=1.4, zorder=2)
            ax.plot(p, yy, "D", color=col, ms=5.0, mfc=col, mec=col, zorder=4)
            ax.text(-0.012, yy, lab, ha="right", va="center", transform=tr,
                    fontsize=7.4, fontweight="bold", clip_on=False)

    ax.axvline(0, color="#404040", lw=0.75, ls=(0, (4, 2.5)), zorder=1)

    vals = []
    for q in ("DI", "SI", "SEC"):
        for r in arms[q]:
            vals += [(math.exp(r["ln"] - 1.959964 * r["se"]) - 1) * 100,
                     (math.exp(r["ln"] + 1.959964 * r["se"]) - 1) * 100]
    lo_x = math.floor(min(vals) / 25) * 25 - 12
    hi_x = math.ceil(max(vals) / 25) * 25 + 12
    ax.set_xlim(lo_x, hi_x)
    ax.set_ylim(-DEPTH - 0.55, 0.80)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.75)
    ax.tick_params(axis="x", width=0.75, length=2.5, labelsize=7.4)
    ax.set_xlabel("Change vs control (%)", fontsize=7.8, labelpad=3)
    ax.set_xticks([t for t in range(int(lo_x // 50 * 50), int(hi_x) + 1, 50)])

    key = [Line2D([0], [0], marker="s", color="#595959", lw=0, ms=3.4,
                  label="single arm (95% CI)"),
           Line2D([0], [0], marker="D", color="#595959", lw=0, ms=4.2,
                  label="pooled (DL, Knapp\u2013Hartung)")]
    ax.legend(handles=key, loc="lower right", frameon=False, fontsize=6.8,
              handletextpad=0.5, labelspacing=0.4, borderpad=0.2)

    fig.subplots_adjust(left=0.335, right=0.985, top=0.985, bottom=0.085)
    save(fig, "Fig2_forest_three_axes_v500")


# ══════════════════════════════════════════════════════════════════
# FIG 3 — 分解（9 臂，Carr 的 DI 为其分项乘积，故仅作闭合点展示）
# ══════════════════════════════════════════════════════════════════
def fig3():
    arms, _ = load()
    si = {(r["study"], r["arm"]): r for r in arms["SI"]}
    sec = {(r["study"], r["arm"]): r for r in arms["SEC"]}
    di = {(r["study"], r["arm"]): r for r in arms["DI"]}

    labels, sv, cv, dv = [], [], [], []
    for st in ORDER:
        for arm in ARMORDER[st]:
            r = di.get((st, arm))
            if r is None:
                continue
            s_row = si.get((st, arm))
            c_row = sec.get((st, arm))
            if s_row is None or c_row is None:
                continue
            labels.append(label(st, arm, r["stratum"]))
            sv.append(s_row["ln"])
            cv.append(c_row["ln"])
            dv.append(r["ln"])

    n = len(labels)
    fig, ax = plt.subplots(figsize=(PRINT_W * CM, (n * 0.72 + 3.0) * CM))
    ys = np.arange(n)[::-1]

    for y, s, c, d in zip(ys, sv, cv, dv):
        ax.barh(y, s, height=0.56, color=C_SI, edgecolor="none", zorder=2)
        start = s if c >= 0 else s + c
        if c >= 0:
            ax.barh(y, abs(c), left=start, height=0.56, color=C_SEC,
                    edgecolor="none", zorder=3)
        else:
            ax.barh(y, abs(c), left=start, height=0.56, color="white",
                    edgecolor=C_SEC, linewidth=0.9, zorder=3)
        ax.plot(d, y, "D", color=C_DI, ms=4.4, mfc="white", mew=1.1, zorder=5)

    ax.axvline(0, color="#404040", lw=0.75, zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=7.2)
    ax.set_ylim(-0.75, n - 0.25)
    ax.set_xlabel("Change vs control (\u0394ln unit)", fontsize=7.8, labelpad=3)
    ax.tick_params(axis="x", width=0.75, length=2.5, labelsize=7.4)
    ax.tick_params(axis="y", length=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    handles = [Patch(facecolor=C_SI, label="insulin sensitivity component"),
               Patch(facecolor=C_SEC, label="insulin secretion component"),
               Line2D([0], [0], marker="D", color="none", mfc="white", mec="black",
                      mew=1.1, ms=4.4, label="observed \u0394ln DI")]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.42, 1.005),
              ncol=3, frameon=False, fontsize=7.0, handlelength=1.3,
              columnspacing=1.4, handletextpad=0.5)

    fig.subplots_adjust(left=0.275, right=0.988, top=0.885, bottom=0.135)
    save(fig, "Fig3_decomposition_v500")


if __name__ == "__main__":
    print("Fig 2 (9 arms, three strata):")
    fig2()
    print("Fig 3 (9 arms):")
    fig3()
