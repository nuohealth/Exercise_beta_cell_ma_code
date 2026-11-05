#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FIG 1 — PRISMA 2020 流程图（v3.18 版式重做，只重建 Fig1）

相对 v3.17 的改动（全部针对「左边文字、连接线、色块」三者位置）：
  1. **左侧阶段带**：原来阶段名直接摆在 0.62 cm 处的空白里，与色块无视觉关系。
     现在划出 1.58 cm 宽的独立阶段带，带浅灰底（#F4F4F4），右缘一条细竖线作分隔，
     阶段名旋转 90° 居中于带内 —— 这是 PRISMA 2020 / SCI 流程图的常规版式。
  2. **色块列位**：左列 x=1.78 w=6.90（右缘 8.68），右列 x=10.00 w=5.60。
     横箭头净长 10.00-8.68 = **1.32 cm**（v3.16 为 0.30，v3.17 为 0.70），杆身清晰。
  3. **色块高度节奏**：设 MIN_H=1.30 cm 下限，并把框内文字块**垂直居中**
     （原实现文字一律贴框顶，单子行行显得上重下空）。
  4. **配色分层**：identification / screening / included 三段各一色系，
     左列实色、右列同色系浅色，避免左右同色导致「排除框」与「主流程框」混淆。
  5. **竖箭头**：起点仍取左盒底（v3.17 修复项，保留），并加 shrinkB 防止箭头尖端
     插入下一框边框。

自检：渲染后逐段取包围盒，越界即报溢量（cm）。
v5.0 数值链：853+477 → 去重 868 → 637 → 431 → 113 → 44 → 主检索纳入 7；
补检索（测量锚，D27）871 → 87 有运动臂 → 64 有全文 → 增纳 1；合计合格 8，入合并 7，供臂 5 项 9 臂。
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(BASE, "figures")
VEC = os.path.join(FIG, "vector")
PRINT_W = 16.0

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["savefig.dpi"] = 300
matplotlib.rcParams["figure.dpi"] = 300

CM = 1 / 2.54

BAND_W = 1.58                    # 左侧阶段带宽度
LX, LW = 1.78, 6.90              # 左列
RX, RW = 10.00, 5.60             # 右列
TOP_PAD, BOT_PAD = 0.52, 0.30
GAP = 0.66
MIN_H = 1.30

C_MAIN = {"id": "#DDE8F4", "scr": "#DDE8F4", "inc": "#DDE8F4"}
C_SIDE = {"id": "#FFFFFF", "scr": "#FFFFFF", "inc": "#FFFFFF"}
INK = "#1A1A1A"


def save(fig, name):
    p_png = os.path.join(FIG, name + ".png")
    fig.savefig(p_png, dpi=300, facecolor="white", transparent=False)
    fig.savefig(os.path.join(VEC, name + ".pdf"), facecolor="white")
    plt.close(fig)
    from PIL import Image
    im = Image.open(p_png)
    if im.info.get("dpi", (0, 0))[0] != 300:
        im.convert("RGB").save(p_png, dpi=(300, 300))
    return p_png


def fig1_prisma():
    F_MAIN, F_SUB, F_STAGE = 7.3, 6.8, 7.3
    LH_MAIN, LH_SUB = 0.325, 0.305
    PAD_T, PAD_B = 0.21, 0.19

    # (左标题, 左子行[], 右标题, 右子行[], 段)
    ROWS = [
        ("Records identified from databases",
         ["PubMed n = 853", "Europe PMC n = 477 (15 new)"],
         "Duplicate records removed", ["n = 462 (Europe PMC only)"], "id"),
        ("Records screened on title/abstract", ["n = 868"],
         "Records excluded",
         ["n = 231",
          "off topic 94 \u00b7 animal 64 \u00b7 acute 19",
          "GDM 17 \u00b7 type 1 diabetes 5 \u00b7 animal+acute 4",
          "no abstract 16 \u00b7 triage only 12"], "scr"),
        ("Records assessed by rule-based screening", ["n = 637"],
         "Records excluded",
         ["n = 206",
          "drug in title 90 \u00b7 not a trial 88",
          "duration < 4 weeks 19 \u00b7 no exercise 5",
          "acute 2 \u00b7 no outcome 2"], "scr"),
        ("Records assessed for a secretion outcome", ["n = 431"],
         "Records excluded",
         ["n = 318", "no insulin-secretion outcome reported"], "scr"),
        ("Reports sought for retrieval", ["n = 113"],
         "Reports not retrieved", ["n = 69 (S5)"], "scr"),
        ("Reports assessed for eligibility", ["n = 44"],
         "Records excluded",
         ["n = 37",
          "no exercise 12 \u00b7 diet or drug only 9",
          "no secretion outcome 8 \u00b7 protocol or review 6",
          "algebraic function of fasting covariates 1",
          "no separable exercise arm 1"], "scr"),
        ("Supplementary search (D27)",
         ["measurement-anchored, n = 871 records",
          "87 described an exercise arm",
          "64 with retrievable full text"],
         "Records excluded",
         ["n = 86",
          "no extractable effect estimate",
          "or criteria unmet (D27)"], "scr"),
        ("Studies included in the review",
         ["n = 8",
          "7 entered the synthesis",
          "5 contributed 9 pooled arms",
          "2 synthesised qualitatively only"],
         "Construct exemplar, outside pooling",
         ["n = 1 \u2014 outcomes are closed-form",
          "functions of fasting covariates"], "inc"),
    ]

    def content_h(n_sub):
        return PAD_T + PAD_B + LH_MAIN + n_sub * LH_SUB

    heights = [max(content_h(len(ls)), content_h(len(rs)), MIN_H)
               for _, ls, _, rs, _ in ROWS]
    total_h = TOP_PAD + sum(heights) + GAP * (len(ROWS) - 1) + BOT_PAD

    fig, ax = plt.subplots(figsize=(PRINT_W * CM, total_h * CM))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_xlim(0, PRINT_W)
    ax.set_ylim(0, total_h)
    ax.axis("off")

    y_content_top = total_h - TOP_PAD
    y_content_bot = BOT_PAD

    # --- 左侧阶段带（浅灰底 + 右缘分隔线）---
    ax.add_patch(Rectangle((0.08, y_content_bot - 0.10), BAND_W - 0.16,
                           (y_content_top - y_content_bot) + 0.20,
                           fc="#F4F4F4", ec="none", zorder=0))
    ax.plot([BAND_W, BAND_W], [y_content_bot - 0.10, y_content_top + 0.10],
            color="#BFBFBF", lw=0.7, zorder=1)

    def draw_box(x, y_top, w, h, title, subs, fc, tag, checks):
        ax.add_patch(FancyBboxPatch((x, y_top - h), w, h,
                                    boxstyle="round,pad=0.06,rounding_size=0.09",
                                    fc=fc, ec=INK, lw=0.75, zorder=2))
        ch = content_h(len(subs))
        off = (h - ch) / 2.0                      # 文字块垂直居中
        cx, yy = x + w / 2.0, y_top - off - PAD_T
        ax.text(cx, yy, title, ha="center", va="top", fontsize=F_MAIN,
                color=INK, zorder=3)
        yy -= LH_MAIN
        for sub in subs:
            ax.text(cx, yy, sub, ha="center", va="top", fontsize=F_SUB,
                    color="#333333", zorder=3)
            yy -= LH_SUB
        checks.append((x, w, y_top, h, tag))

    def varrow(x, y1, y2):
        ax.add_patch(FancyArrowPatch((x, y1), (x, y2), arrowstyle="-|>",
                                     mutation_scale=7.5, color="#3A3A3A",
                                     lw=0.85, shrinkA=0, shrinkB=1.2, zorder=1))

    def harrow(x1, y, x2):
        ax.add_patch(FancyArrowPatch((x1, y), (x2, y), arrowstyle="-|>",
                                     mutation_scale=7.5, color="#3A3A3A",
                                     lw=0.85, shrinkA=0, shrinkB=1.2, zorder=1))

    # 行带高度 = max(左内容, 右内容, MIN_H)；左框按内容定高并居中于行带
    checks, spans, geom = [], [], []
    y_top = y_content_top
    for idx, (lt, ls, rt, rs, seg) in enumerate(ROWS):
        rh = heights[idx]
        hl = max(content_h(len(ls)), MIN_H)
        hr = max(content_h(len(rs)), MIN_H) if rt else 0.0
        off_l = (rh - hl) / 2.0
        off_r = (rh - hr) / 2.0
        spans.append((y_top, rh, off_l, hl))
        geom.append((y_top - off_l, hl))
        draw_box(LX, y_top - off_l, LW, hl, lt, ls, C_MAIN[seg],
                 f"\u884c{idx+1}\u5de6", checks)
        if rt:
            draw_box(RX, y_top - off_r, RW, hr, rt, rs, C_SIDE[seg],
                     f"\u884c{idx+1}\u53f3", checks)
            harrow(LX + LW, y_top - off_l - hl / 2.0, RX)
        y_top -= (rh + GAP)

    # --- 阶段名 ---
    # 竖向箭头：连接相邻左框（框底 → 下一框顶）
    for i in range(len(ROWS) - 1):
        y1 = geom[i][0] - geom[i][1]
        y2 = geom[i + 1][0]
        varrow(LX + LW / 2.0, y1, y2)

    STAGES = [(0, 0, "Identification"), (1, 6, "Screening"), (7, 7, "Included")]
    for i0, i1, label in STAGES:
        y_hi = spans[i0][0]
        y_lo = spans[i1][0] - spans[i1][1]
        ax.text(BAND_W / 2.0, (y_hi + y_lo) / 2.0, label, rotation=90,
                ha="center", va="center", fontsize=F_STAGE, fontweight="bold",
                color=INK, zorder=3)

    # --- 自检：文字必须落在所属框内 ---
    fig.canvas.draw()
    inv = ax.transData.inverted()
    bad = []
    for (bx, bw, by, bh, tag) in checks:
        for t in ax.texts:
            tx, ty = t.get_position()
            if not (bx - 0.01 <= tx <= bx + bw + 0.01 and by - bh - 0.01 <= ty <= by + 0.01):
                continue
            bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
            x0, y0 = inv.transform((bb.x0, bb.y0))
            x1, y1 = inv.transform((bb.x1, bb.y1))
            if x0 < bx - 0.02 or x1 > bx + bw + 0.02:
                bad.append(f"{tag} | {t.get_text()[:44]!r} \u5de6{ x0-bx:+.2f} \u53f3{x1-(bx+bw):+.2f} cm")
            if y1 > by + 0.02 or y0 < by - bh - 0.02:
                bad.append(f"{tag} | {t.get_text()[:44]!r} \u4e0a{y1-by:+.2f} \u4e0b{y0-(by-bh):+.2f} cm")
    print("  [自检] 文字越界：" if bad else "  [自检] 全部文字在框内（含上下缘）")
    for b in bad:
        print("    " + b)

    p = save(fig, "Fig1_PRISMA_v500")
    print(f"Fig1 重建: {p}")
    print(f"  画布 {PRINT_W:.1f} x {total_h:.2f} cm | 阶段带 {BAND_W} cm | "
          f"横箭头净长 {RX-(LX+LW):.2f} cm")
    for i, h in enumerate(heights):
        print(f"  行{i+1} 高 {h:.2f} cm")


if __name__ == "__main__":
    fig1_prisma()
