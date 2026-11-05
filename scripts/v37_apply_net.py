#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v3.7 步1 —— 把 Slentz 2009 三个运动臂统一到 net（相对对照）口径

背景
----
v3.6 稿件内部混用两种尺度：
  · DOSE-EX 三臂 / Karstoft —— 本来就是组间对比（net）
  · Slentz 三臂（Table 2 + Fig 3）—— 用的是 from-baseline
  · 但 Table 3 合并、Table 2 脚注("net of the control arm")、Fig 2 横轴("vs control")
    都已是 net 口径
→ 同一条臂在图里一个数、在合并里另一个数，审稿人一验即穿。

本脚本把 Slentz 三臂改为 net = ln(arm_post/arm_base) − ln(CON_post/CON_base)。
标准误按独立近似合成：se_net = sqrt(se_arm^2 + se_con^2)（VIF 约 2，
即每臂方差膨胀 ~2 倍，见 p10_variance_inflation.py）。

CON 臂本身不是效应量（它是参照），保留在 CSV 里但 include_in_fig=0，
其 8 个月自身漂移写入脚注，供审稿人核对。

单一数据源：analysis/decomposition.csv
输出：原地更新 + analysis/decomposition_v36_backup.csv（首次运行留档）
"""
import csv, math, os, shutil

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P = os.path.join(BASE, "analysis", "decomposition.csv")
BK = os.path.join(BASE, "analysis", "decomposition_v36_backup.csv")

rows = list(csv.DictReader(open(P)))
cols = list(rows[0].keys())

# 首次运行留档（可复现性：改前状态必须可追）
if not os.path.exists(BK):
    shutil.copy(P, BK)

def f(x):
    return float(x) if x not in ("", None) else None

con = next(r for r in rows if r["study"] == "Slentz 2009" and r["arm"] == "CON")
con_di, con_si, con_se = f(con["dln_DI"]), f(con["dln_Si"]), f(con["dln_Sec"])
con_sedi, con_sesi, con_sesec = f(con["se_DI"]), f(con["se_Si"]), f(con["se_Sec"])

changed = []
for r in rows:
    if r["study"] != "Slentz 2009":
        continue
    if r["arm"] == "CON":
        r["scale"] = "reference (drift from baseline, 8 mo)"
        r["include_in_fig"] = "0"
        r["note"] = (f'CON 自身漂移（非效应量）: DI {con_di:+.4f} / Si {con_si:+.4f} / '
                     f'Sec {con_se:+.4f} ln；AIRg 491→513, Si 3.5→3.1, DI 1562→1497')
        continue

    for q, cq, csq in (("DI", con_di, con_sedi), ("Si", con_si, con_sesi), ("Sec", con_se, con_sesec)):
        a, sa = f(r["dln_" + q]), f(r["se_" + q])
        if a is None:
            continue
        r["dln_" + q] = round(a - cq, 6)
        r["se_" + q] = round(math.sqrt(sa ** 2 + csq ** 2), 6)

    di, si, sec = f(r["dln_DI"]), f(r["dln_Si"]), f(r["dln_Sec"])
    resid = di - si - sec
    r["residual"] = round(resid, 6)
    r["resid_se"] = round(math.sqrt(f(r["se_DI"]) ** 2 + f(r["se_Si"]) ** 2 + f(r["se_Sec"]) ** 2), 6)
    r["share_Si"] = round(si / di, 4) if di else ""
    r["share_Sec"] = round(sec / di, 4) if di else ""
    r["scale"] = "arithmetic (net of control)"
    r["include_in_fig"] = "1"
    r["note"] = (r["note"].split(" | AIRg")[0]
                 + f" | net of CON: DI {(math.exp(di)-1)*100:+.1f}% "
                   f"Si {(math.exp(si)-1)*100:+.1f}% Sec {(math.exp(sec)-1)*100:+.1f}%")
    changed.append((r["arm"], di, si, sec, resid))

# 非 Slentz 行补列
for r in rows:
    if "include_in_fig" not in r:
        r["include_in_fig"] = "1"

cols = list(dict.fromkeys(cols + ["include_in_fig"]))
with open(P, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    for r in rows:
        w.writerow({c: r.get(c, "") for c in cols})

print("Slentz 三臂 → net 完成")
print(f"{'arm':10s} {'dln_DI':>9s} {'dln_Si':>9s} {'dln_Sec':>9s} {'resid':>8s}  (旧残差)")
for arm, di, si, sec, resid in changed:
    old = {"LOWMOD": 0.0450, "LOWVIG": 0.0596, "HIGHVIG": 0.0482}[arm]
    print(f"{arm:10s} {di:9.4f} {si:9.4f} {sec:9.4f} {resid:8.4f}  (was {old:+.4f})")
print(f"\nCON 漂移（参照，不入图）: DI {con_di:+.4f} Si {con_si:+.4f} Sec {con_se:+.4f}")
