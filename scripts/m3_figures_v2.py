#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Diabetologia 主图重制 v2 —— 按终稿印刷尺寸设计画布（16.0 cm 版心宽）

设计原则
  1. 画布宽度 = 终稿印刷宽度（16.0 cm），插入 docx 时 1:1 不缩放
     → 脚本里写的字号 = 印出来的字号，不再被二次压缩
  2. 无图内标题、无图内统计值框（期刊要求：omit headings, captions and statistical values）
  3. 配色明度分离：敏感性深蓝 / 分泌橙，灰度对比度 ≈3:1，黑白打印可辨
  4. 线宽 0.75 pt（期刊推荐），轴黑色，无网格、无边框
  5. 字体 Liberation Sans（与 Arial 度量兼容，期刊要求 sans serif）

输出：figures/*.png (300 dpi, 白底 RGB) + figures/vector/*.pdf (矢量)
"""
import json, math, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.transforms import blended_transform_factory
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(BASE, "figures")
VEC = os.path.join(FIG, "vector")
os.makedirs(FIG, exist_ok=True)
os.makedirs(VEC, exist_ok=True)

CM = 1 / 2.54                      # cm -> inch
PRINT_W = 16.0                     # 版心宽（A4 21cm − 2×2.5cm 页边距）

plt.rcParams.update({
    "font.family": "Liberation Sans",
    "font.size": 8,
    "axes.linewidth": 0.75,
    "axes.edgecolor": "black",
    "savefig.dpi": 300,
    "figure.dpi": 300,
    "savefig.facecolor": "white",
    "savefig.transparent": False,
    "pdf.fonttype": 42,            # TrueType 嵌入，矢量 PDF 文字可编辑
    "ps.fonttype": 42,
})

# 明度分离配色（灰度可辨）
C_SI  = "#1F4E79"   # 深蓝 —— 胰岛素敏感性
C_SEC = "#E08214"   # 橙   —— 胰岛素分泌
C_DI  = "#000000"   # 黑   —— 处置指数
C_GREY = "#8C8C8C"
C_BOX = "#F2F2F2"
C_BOX_EX = "#FBE9D7"


def save(fig, name, vector=True):
    """不裁切：画布宽度严格等于 16.0 cm，插入 docx 时 1:1 不缩放。"""
    from PIL import Image
    p_png = os.path.join(FIG, name + ".png")
    fig.savefig(p_png, dpi=300, facecolor="white", transparent=False)
    if vector:
        fig.savefig(os.path.join(VEC, name + ".pdf"), facecolor="white")
    plt.close(fig)

    # 印刷需白底 RGB，去掉 alpha 通道
    im = Image.open(p_png)
    if im.mode != "RGB":
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[-1] if im.mode == "RGBA" else None)
        bg.save(p_png, dpi=(300, 300))
    print(f"  wrote {name}.png" + (" + vector pdf" if vector else ""))


def clean_axes(ax, spines=("top", "right")):
    for s in spines:
        ax.spines[s].set_visible(False)
    ax.tick_params(width=0.75, length=2.5)


# ══════════════════════════════════════════════════════════════════
# FIG 1 — PRISMA 2020 flow diagram
# 每个节点数字已与正文/补充材料勾稽：868 → 637 → 431 → 113 → 44 → 6
# ══════════════════════════════════════════════════════════════════
def fig1_prisma():
    """
    PRISMA 2020 流程图。
    布局引擎：按「行」堆叠，行高 = max(左盒, 右盒)，行间留固定空隙
    → 结构上不可能出现盒/文字重叠（不靠手算坐标）。
    数值链已与正文、ESM 勾稽：
      853+477 → 去重 868 → 筛 637 → 规则筛 431 → 核心 113 → 全文 44 → 纳入 6
    """
    F_MAIN, F_SUB, F_SML = 7.2, 6.7, 6.5
    LH_MAIN = 0.32          # cm per line
    LH_SUB = 0.305
    PAD_T, PAD_B = 0.22, 0.20
    GAP = 0.62              # 行间空隙（箭头空间）

    LX, LW = 2.4, 7.3       # 左列
    RX, RW = 10.0, 5.6      # 右列（右缘留 0.4 cm 边距，防边框被画布裁切）

    # (左标题, 左子行[], 右标题, 右子行[]) —— 任一可为 None
    ROWS = [
        ("Records identified from databases",
         ["PubMed n = 853", "Europe PMC n = 477 (15 new)"],
         "Duplicate records removed", ["n = 462 (Europe PMC only)"]),

        ("Records screened on title/abstract",
         ["n = 868"],
         "Records excluded",
         ["n = 231",
          "off topic 94 \u00b7 animal 64 \u00b7 acute 19",
          "GDM 17 \u00b7 type 1 diabetes 5 \u00b7 animal+acute 4",
          "no abstract 16 \u00b7 triage only 12"]),

        ("Records assessed by rule-based screening",
         ["n = 637"],
         "Records excluded",
         ["n = 206",
          "drug in title 90 \u00b7 not a trial 88",
          "duration < 4 weeks 19 \u00b7 no exercise 5",
          "acute 2 \u00b7 no outcome 2"]),

        ("Records assessed for a secretion outcome",
         ["n = 431"],
         "Records excluded",
         ["n = 318", "no insulin-secretion outcome reported"]),

        ("Reports sought for retrieval",
         ["n = 113"],
         "Reports not retrieved",
         ["n = 69 (ESM Table S5)"]),

        ("Reports assessed for eligibility",
         ["n = 44"],
         "Reports excluded",
         ["n = 37",
          "no exercise 12 \u00b7 diet or drug only 9",
          "no secretion outcome 8 \u00b7 protocol or review 6",
          "algebraic function of fasting covariates 1",
          "no separable exercise arm 1"]),

        ("Studies included in review",
         ["n = 6",
          "4 contributed 8 arms to pooling",
          "2 synthesised qualitatively"],
         "Construct exemplar, outside pooling",
         ["n = 1 \u2014 outcomes are closed-form",
          "functions of fasting covariates"]),
    ]

    def box_h(title, subs):
        h = PAD_T + PAD_B + LH_MAIN
        if subs:
            h += len(subs) * LH_SUB
        return h

    # 计算每行高度与总高
    heights = []
    for lt, ls, rt, rs in ROWS:
        hl = box_h(lt, ls) if lt else 0
        hr = box_h(rt, rs) if rt else 0
        heights.append(max(hl, hr))

    total_h = sum(heights) + GAP * (len(ROWS) - 1) + 0.5

    W, H = PRINT_W * CM, total_h * CM
    fig, ax = plt.subplots(figsize=(W, H))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_xlim(0, 16.0)
    ax.set_ylim(0, total_h)
    ax.axis("off")

    def draw_box(x, y_top, w, title, subs, fc):
        h = box_h(title, subs)
        ax.add_patch(FancyBboxPatch((x, y_top - h), w, h,
                                    boxstyle="round,pad=0.10,rounding_size=0.10",
                                    fc=fc, ec="black", lw=0.7, zorder=2))
        cx = x + w / 2.0
        yy = y_top - PAD_T
        ax.text(cx, yy, title, ha="center", va="top", fontsize=F_MAIN,
                zorder=3)
        yy -= LH_MAIN
        for sub in subs:
            ax.text(cx, yy, sub, ha="center", va="top", fontsize=F_SUB,
                    color="#222222", zorder=3)
            yy -= LH_SUB
        return h

    def varrow(x, y1, y2):
        ax.add_patch(FancyArrowPatch((x, y1), (x, y2), arrowstyle="-|>",
                                     mutation_scale=7, color="#404040",
                                     lw=0.8, shrinkA=0, shrinkB=0, zorder=1))

    def harrow(x1, y, x2):
        ax.add_patch(FancyArrowPatch((x1, y), (x2, y), arrowstyle="-|>",
                                     mutation_scale=7, color="#404040",
                                     lw=0.8, shrinkA=0, shrinkB=0, zorder=1))

    y_top = total_h - 0.25
    spans = []                                  # (y_top, row_h) per row, for stage labels
    for idx, (lt, ls, rt, rs) in enumerate(ROWS):
        row_h = heights[idx]                    # 行高 = max(左,右)
        spans.append((y_top, row_h))
        hl = draw_box(LX, y_top, LW, lt, ls,
                      "#D8EBD5" if idx >= 4 else "#DCE9F5")
        if rt:
            draw_box(RX, y_top, RW, rt, rs,
                     "#EFE3D2" if idx == 6 else "#FBE9D7")
            harrow(LX + LW, y_top - hl / 2.0, RX)

        if idx < len(ROWS) - 1:
            varrow(LX + LW / 2.0, y_top - row_h, y_top - row_h - GAP)
            y_top = y_top - row_h - GAP
        else:
            y_top = y_top - row_h

    # PRISMA 2020 左侧阶段标签（Identification / Screening / Included）
    STAGES = [(0, 0, "Identification"), (1, 5, "Screening"), (6, 6, "Included")]
    for i0, i1, label in STAGES:
        y_hi = spans[i0][0]
        y_lo = spans[i1][0] - spans[i1][1]
        ax.text(0.62, (y_hi + y_lo) / 2.0, label, rotation=90,
                ha="center", va="center", fontsize=F_MAIN, fontweight="bold",
                color="#333333")

    save(fig, "Fig1_PRISMA")


# ══════════════════════════════════════════════════════════════════
# FIG 2 — 三轴森林图（单一坐标，同尺度并列）
# 关键读法：DI 与敏感性右移，分泌停在零线 —— 一眼可见
# ══════════════════════════════════════════════════════════════════
def fig2_forest():
    A = json.load(open(os.path.join(BASE, "analysis", "meta_A.json")))
    per = A["per_arm_data"]
    pooled = A["A_exercise_vs_control"]

    def pct(ln):
        return (math.exp(ln) - 1) * 100

    # 面板顺序 = 论文论证顺序 DI → SI → SEC，与图注 a/b/c 一一对应
    sections = [
        ("DI",  "Disposition index",   C_DI,  "a"),
        ("SI",  "Insulin sensitivity", C_SI,  "b"),
        ("SEC", "Insulin secretion",   C_SEC, "c"),
    ]

    ROW = 0.33          # 每行 0.33 cm（印刷高度，非画布缩放）
    HDR_GAP = 0.85
    POOL_GAP = 1.05

    rows = []      # (y, kind, payload)
    y = 0.0
    for q, title, col, pletter in sections:
        rows.append((y, "header", (title, col, pletter)))
        y -= HDR_GAP
        arms = [r for r in per if r["q"] == q]
        for r in arms:
            lo = pct(r["ln"] - 1.959964 * r["se"])
            hi = pct(r["ln"] + 1.959964 * r["se"])
            rows.append((y, "arm", (f'{r["study"]}  ·  {r["arm"]}', pct(r["ln"]), lo, hi, col)))
            y -= ROW
        m = pooled[q]
        rows.append((y, "pooled", ("Pooled (k=%d)" % m["k"], m["pct"],
                                   m["ci_lo_pct"], m["ci_hi_pct"], col)))
        y -= POOL_GAP

    DEPTH = -y
    fig_h = DEPTH + 1.1
    fig, ax = plt.subplots(figsize=(PRINT_W * CM, fig_h * CM))

    tr = blended_transform_factory(ax.transAxes, ax.transData)

    for (yy, kind, payload) in rows:
        if kind == "header":
            title, col, pletter = payload
            ax.text(-0.012, yy, title, ha="right", va="center",
                    transform=tr, fontsize=8.4, fontweight="bold", color=col,
                    clip_on=False)
            # 分图字母（期刊体例：粗体，置于面板左上）
            ax.text(-0.012, yy + 0.34, pletter, ha="right", va="center",
                    transform=tr, fontsize=8.4, fontweight="bold",
                    color="black", clip_on=False)
            ax.axhline(yy - 0.50, color=col, lw=0.75, alpha=0.45,
                       xmin=0, xmax=1)
        elif kind == "arm":
            lab, p, lo, hi, col = payload
            ax.plot([lo, hi], [yy, yy], color=col, lw=1.0,
                    solid_capstyle="butt", zorder=2)
            ax.plot(p, yy, "s", color=col, ms=3.6, zorder=3)
            ax.text(-0.012, yy, lab, ha="right", va="center",
                    transform=tr, fontsize=7.2, clip_on=False)
        else:
            lab, p, lo, hi, col = payload
            ax.plot([lo, hi], [yy, yy], color=col, lw=1.4, zorder=2)
            ax.plot(p, yy, "D", color=col, ms=5.0, mfc=col, mec=col, zorder=4)
            ax.text(-0.012, yy, lab, ha="right", va="center", transform=tr,
                    fontsize=7.4, fontweight="bold", clip_on=False)

    ax.axvline(0, color="#404040", lw=0.75, ls=(0, (4, 2.5)), zorder=1)

    ax.set_xlim(-44, 238)
    ax.set_ylim(-DEPTH - 0.55, 0.80)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.75)
    ax.tick_params(axis="x", width=0.75, length=2.5, labelsize=7.4)
    ax.set_xlabel("Change vs control (%)", fontsize=7.8, labelpad=3)
    # 零线与刻度
    ax.set_xticks([-25, 0, 50, 100, 150, 200])


    # 符号说明：单臂=方块，合并=菱形
    from matplotlib.lines import Line2D
    key = [
        Line2D([0], [0], marker="s", color="#595959", lw=0, ms=3.4,
               label="single arm (95% CI)"),
        Line2D([0], [0], marker="D", color="#595959", lw=0, ms=4.2,
               label="pooled (DL, Knapp–Hartung)"),
    ]
    ax.legend(handles=key, loc="lower right", frameon=False, fontsize=6.8,
              handletextpad=0.5, labelspacing=0.4, borderpad=0.2)

    fig.subplots_adjust(left=0.315, right=0.985, top=0.985, bottom=0.085)
    save(fig, "Fig2_forest_three_axes")


# ══════════════════════════════════════════════════════════════════
# FIG 3 — 分解：每个臂的 DI 变化拆成敏感性 + 分泌
# 条从零起算，敏感性段末端接分泌段，菱形 = 实测 DI
# ══════════════════════════════════════════════════════════════════
def fig3_decomp():
    import csv
    D = [r for r in csv.DictReader(open(os.path.join(BASE, "analysis", "decomposition.csv")))
         if r["dln_DI"] not in ("", None) and r.get("include_in_fig", "1") == "1"]

    order = {"Slentz 2009": 0, "DOSE-EX 2023": 1, "Karstoft 2014": 2}
    arm_order = {"LOWMOD": 0, "LOWVIG": 1, "HIGHVIG": 2,
                 "DCON": 0, "MED": 1, "HED": 2, "IWT vs CON": 0}
    D.sort(key=lambda r: (order.get(r["study"], 9), arm_order.get(r["arm"], 9)))

    labels, si, sec, di = [], [], [], []
    for r in D:
        lab = r["arm"].replace(" vs CON", "")
        labels.append(f'{r["study"]}  ·  {lab}')
        si.append(float(r["dln_Si"]))
        sec.append(float(r["dln_Sec"]))
        di.append(float(r["dln_DI"]))

    n = len(labels)
    fig, ax = plt.subplots(figsize=(PRINT_W * CM, (n * 0.72 + 3.0) * CM))
    ys = np.arange(n)[::-1]

    for y, s, c, d in zip(ys, si, sec, di):
        # 敏感性段：0 → s
        ax.barh(y, s, height=0.56, color=C_SI, edgecolor="none", zorder=2)
        # 分泌段：从敏感性命中末端延伸（正负皆按几何拼接）
        start = s if c >= 0 else s + c
        # 负分泌段：白色填充 + 橙色描边，从敏感性末端向左"咬掉"一块
        # → 语义正确（分泌为负 = 从 DI 中扣除），且不会被蓝条吞掉
        if c >= 0:
            ax.barh(y, abs(c), left=start, height=0.56, color=C_SEC,
                    edgecolor="none", zorder=3)
        else:
            ax.barh(y, abs(c), left=start, height=0.56, color="white",
                    edgecolor=C_SEC, linewidth=0.9, zorder=3)
        # 实测 DI（菱形）
        ax.plot(d, y, "D", color=C_DI, ms=4.4, mfc="white", mew=1.1, zorder=5)

    ax.axvline(0, color="#404040", lw=0.75, zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=7.2)
    ax.set_ylim(-0.75, n - 0.25)
    ax.set_xlabel("Change vs control (Δln unit)", fontsize=7.8, labelpad=3)
    ax.tick_params(axis="x", width=0.75, length=2.5, labelsize=7.4)
    ax.tick_params(axis="y", length=0)
    clean_axes(ax, ("top", "right"))

    # 图例（无框、置于底部）
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor=C_SI, label="insulin sensitivity component"),
        Patch(facecolor=C_SEC, label="insulin secretion component"),
        Line2D([0], [0], marker="D", color="none", mfc="white", mec="black",
               mew=1.1, ms=4.4, label="observed Δln DI"),
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.42, 1.005),
              ncol=3, frameon=False, fontsize=7.0, handlelength=1.3,
              columnspacing=1.4, handletextpad=0.5)

    fig.subplots_adjust(left=0.255, right=0.988, top=0.885, bottom=0.135)
    save(fig, "Fig3_decomposition")


# ══════════════════════════════════════════════════════════════════
# FIG 4 — DOSE-EX：a) 剂量反应  b) 等减重对比（决定性证据）
# ══════════════════════════════════════════════════════════════════
def fig4_doseex():
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(PRINT_W * CM, 8.8 * CM),
        gridspec_kw={"width_ratios": [1.28, 1.0], "wspace": 0.34})

    # ── a) 相对 CON 的剂量反应 ──
    arms = ["CON", "DCON", "MED", "HED"]
    DI  = [0, 58, 105, 137]
    ISI = [0, 18, 50, 83]
    ISR = [0, 33, 38, 28]
    xs = np.arange(len(arms))
    # CON（标准照护）只作参照点：空心符号、不连线。
    # 真正的运动剂量轴只有 DCON(0) → MED(3×/wk) → HED(6×/wk)，
    # 把 CON 连进去等于把限食效应混进"剂量"，是逻辑错误。
    ax1.plot([0], [0], "o", color=C_SEC, ms=4.2, mfc="white", mew=1.1, zorder=3)
    ax1.plot([0], [0], "s", color=C_SI,  ms=4.0, mfc="white", mew=1.1, zorder=3)
    ax1.plot([0], [0], "^", color=C_DI,  ms=4.4, mfc="white", mew=1.1, zorder=3)
    dose_x = [1, 2, 3]
    ax1.plot(dose_x, ISR[1:], "-o", color=C_SEC, lw=1.3, ms=4.2, label="Insulin secretion (ISR)")
    ax1.plot(dose_x, ISI[1:], "-s", color=C_SI, lw=1.3, ms=4.0, label="Insulin sensitivity (ISI)")
    ax1.plot(dose_x, DI[1:],  "-^", color=C_DI, lw=1.3, ms=4.4, label="Disposition index (DI)")
    # 剂量轴的虚线引导：CON 与 DCON 之间断开，表示参照点→剂量轴起点的跳跃
    ax1.axvline(0.5, color="#BFBFBF", lw=0.6, ls=(0, (2, 3)), zorder=1)
    ax1.axhline(0, color="#404040", lw=0.75, ls=(0, (4, 2.5)))
    ax1.set_xticks(xs)
    ax1.set_xticklabels(["CON\n(reference)", "DCON\n(0×/wk)", "MED\n(3×/wk)", "HED\n(6×/wk)"],
                        fontsize=7.0)
    ax1.set_xlim(-0.35, 3.35)
    ax1.set_ylim(-16, 152)
    ax1.set_ylabel("Change vs CON (%)", fontsize=7.8, labelpad=3)
    ax1.set_xlabel("Exercise dose (sessions per week, added to caloric restriction)",
                   fontsize=7.6, labelpad=3)
    ax1.text(0.02, 0.96, "a", transform=ax1.transAxes, fontsize=10,
             fontweight="bold", va="top")
    ax1.tick_params(axis="both", width=0.75, length=2.5, labelsize=7.4)
    clean_axes(ax1, ("top", "right"))
    ax1.legend(loc="upper left", bbox_to_anchor=(0.02, 0.90), frameon=False,
               fontsize=6.9, handlelength=1.5, labelspacing=0.35)

    # ── b) 等减重：HED / MED vs DCON（带 95% CI）──
    contrasts = [
        ("DI",  [(50, 10, 104), (29, -5, 77)],  C_DI),
        ("ISI", [(55, 15, 109), (27, -7, 71)],  C_SI),
        ("ISR", [(-4, -15, 9),  (4, -8, 18)],   C_SEC),
    ]
    ypos = {"DI": 2.0, "ISI": 1.0, "ISR": 0.0}
    off = 0.17
    for lab, bars, col in contrasts:
        for i, (p, lo, hi) in enumerate(bars):
            yy = ypos[lab] + (off if i == 0 else -off)
            ax2.plot([lo, hi], [yy, yy], color=col, lw=1.1, solid_capstyle="butt", zorder=2)
            ax2.plot(p, yy, "o" if i == 0 else "o", color=col,
                     ms=4.0 if i == 0 else 3.4, mfc=col if i == 0 else "white",
                     mew=0.0 if i == 0 else 1.0, zorder=3)
    ax2.axvline(0, color="#404040", lw=0.75, ls=(0, (4, 2.5)), zorder=1)
    ax2.set_yticks([2.0, 1.0, 0.0])
    ax2.set_yticklabels(["DI", "ISI", "ISR"], fontsize=7.6)
    ax2.set_ylim(-0.62, 2.62)
    ax2.set_xlim(-30, 125)
    ax2.set_xlabel("Change vs DCON (%)", fontsize=7.8, labelpad=3)
    ax2.tick_params(axis="both", width=0.75, length=2.5, labelsize=7.4)
    ax2.tick_params(axis="y", length=0)
    clean_axes(ax2, ("top", "right"))
    ax2.text(0.02, 0.96, "b", transform=ax2.transAxes, fontsize=10,
             fontweight="bold", va="top")

    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color=C_GREY, mfc=C_GREY, mew=0, ms=4.0,
               label="high dose (6×/wk)"),
        Line2D([0], [0], marker="o", color=C_GREY, mfc="white", mew=1.0, ms=3.4,
               label="moderate dose (3×/wk)"),
    ]
    ax2.legend(handles=handles, loc="lower right", frameon=False, fontsize=6.8,
               handlelength=1.0, labelspacing=0.35)

    fig.subplots_adjust(left=0.075, right=0.985, top=0.955, bottom=0.135)
    save(fig, "Fig4_DOSEEX_armwise")


if __name__ == "__main__":
    print("重制主图（画布 = 印刷尺寸 16.0 cm）")
    fig1_prisma()
    fig2_forest()
    fig3_decomp()
    fig4_doseex()
    print("\n完成。校验：")
    for n in ("Fig1_PRISMA", "Fig2_forest_three_axes", "Fig3_decomposition",
              "Fig4_DOSEEX_armwise"):
        from PIL import Image
        im = Image.open(os.path.join(FIG, n + ".png"))
        print(f"  {n}.png  {im.size[0]}x{im.size[1]} px @300dpi "
              f"= {im.size[0]/300*2.54:.2f} x {im.size[1]/300*2.54:.2f} cm  mode={im.mode}")
