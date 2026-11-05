#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
m2_stratified.py — 按基线血糖状态分层 + Lee 2012 纳入（v3.13 方案 B）

设计依据
--------
v3.12 的合并只用 6 臂 / 3 研究，但其中 Slentz 2009 的三臂来自"肥胖+血脂异常、
糖代谢正常"人群（原文排除糖尿病，基线 FPG 91-94 mg/dl，自述 normal glucose tertiles），
与 §2.3「obesity **with dysglycaemia**」不符；同为正常血糖的 HEARTY 已被排除（S4 #13）。
只剔 Lee 会留下 Slentz，构成同一标准的双标。故改为**分层全报**：

  L1  T2D 层            : DOSE-EX HED/MED, Karstoft IWT
  L2  正常血糖肥胖层    : Slentz LOWMOD/LOWVIG/HIGHVIG, Lee AE/RE
  ALL 全域（两层合并）

Lee 2012 的数只是**相对对照的组间差值（net）**，非组内变化：
  Table 3 脚注原文 "Change values for the intervention groups are the difference
  as compared with control with SE."；对照组的 change 为 imputed means。
  故 ln 尺度用 ln(1 + Δ/B_arm)，SE 由 delta method: SE_Δ / (B_arm + Δ)。
  这与 Slentz 的 (臂内 ln 比值 − CON 臂内 ln 比值) 同为 net 尺度，但实现不同，
  差异已在 S7 敏感性分析中披露。

单一数据源：extraction/pdf/lee_tables/tbl1.jpg(基线) + tbl3.jpg(Δ vs control)
            + scripts/m2_meta.py 内置 DOSE-EX / Karstoft / Slentz 常数

用法: python3 scripts/m2_stratified.py
输出: analysis/stratified.json, analysis/stratified.md
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m2_meta as M

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ---- Lee 2012 Table 1 基线（mean ± SD）与 Table 3 ITT Δ vs control（mean ± SE）----
# 行: (基线 control, 基线 AE, 基线 RE);  Δ: (AE Δ, AE SE, RE Δ, RE SE)
LEE = {
    "DI": {"base": {"AE": 538.9, "RE": 545.5},
           "delta": {"AE": (105.7, 78.0), "RE": (37.4, 77.1)}},
    "SI": {"base": {"AE": 2.2, "RE": 2.9},
           "delta": {"AE": (0.4, 0.2), "RE": (0.8, 0.2)}},
    # first-phase insulin：钳夹动态分泌指标，与 Slentz 的 AIRg 同类
    "SEC": {"base": {"AE": 292.7, "RE": 213.9},
            "delta": {"AE": (3.9, 22.4), "RE": (-38.2, 21.5)}},
}

T2D_STUDIES = {"DOSE-EX 2023", "Karstoft 2014"}
NORM_STUDIES = {"Slentz 2009"}


def lee_arms():
    out = []
    for q, spec in LEE.items():
        for arm in ["AE", "RE"]:
            b = spec["base"][arm]
            d, se_d = spec["delta"][arm]
            post = b + d
            if b <= 0 or post <= 0:
                continue
            out.append({"study": "Lee 2012", "arm": arm, "q": q,
                        "ln": math.log(post / b), "se": se_d / post,
                        "src": "Table 3 (Δ vs control, ITT)"})
    return out


def layer(records, q, studies):
    return [r for r in records if r["q"] == q and r["study"] in studies]


def run():
    base = M.build_set_A()
    allrec = base + lee_arms()

    layers = {
        "L1_T2D": T2D_STUDIES,
        "L2_normoglycaemic_obesity": NORM_STUDIES | {"Lee 2012"},
        "ALL": T2D_STUDIES | NORM_STUDIES | {"Lee 2012"},
    }

    res = {}
    for lname, studies in layers.items():
        res[lname] = {}
        for q in ["DI", "SI", "SEC"]:
            rows = layer(allrec, q, studies)
            if len(rows) < 2:
                res[lname][q] = {"k": len(rows), "pooled": False,
                                 "arms": [(r["study"] + "/" + r["arm"]) for r in rows]}
                continue
            g = M.dl_kh([r["ln"] for r in rows], [r["se"] for r in rows])
            res[lname][q] = {
                "k": g["k"], "pooled": True,
                "pct": round(g["pct"], 1),
                "ci": [round(g["ci_lo_pct"], 1), round(g["ci_hi_pct"], 1)],
                "p": float(f"{g['p']:.4g}"),
                "I2": round(g["I2"], 1),
                "tau2": round(g["tau2"], 5),
                "arms": [(r["study"] + "/" + r["arm"]) for r in rows],
            }
    return res, allrec


def md_table(res):
    order = [("L1_T2D", "T2D"), ("L2_normoglycaemic_obesity", "Normoglycaemic obesity"), ("ALL", "All")]
    lines = ["| Layer | DI | SI | SEC |", "|---|---|---|---|"]
    for key, label in order:
        cells = []
        for q in ["DI", "SI", "SEC"]:
            d = res[key][q]
            if not d["pooled"]:
                cells.append(f"k={d['k']}, not pooled")
            else:
                cells.append(f"**{d['pct']:+.1f}%** ({d['ci'][0]:+.1f}, {d['ci'][1]:+.1f}); "
                             f"p={d['p']}; I²={d['I2']}%; k={d['k']}")
        lines.append(f"| {label} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


if __name__ == "__main__":
    res, allrec = run()
    os.makedirs(os.path.join(ROOT, "analysis"), exist_ok=True)
    with open(os.path.join(ROOT, "analysis/stratified.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)

    tbl = md_table(res)
    with open(os.path.join(ROOT, "analysis/stratified.md"), "w", encoding="utf-8") as f:
        f.write("# Stratified pooled estimates (v3.13)\n\n" + tbl + "\n")

    print(tbl)
    print("\n--- arms per layer ---")
    for k, v in res.items():
        for q in ["DI", "SI", "SEC"]:
            print(f"{k:26s} {q:4s} k={v[q]['k']}  {', '.join(v[q]['arms'])}")
