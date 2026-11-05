#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
m2_v313.py — v3.13（方案 B：分层全报）的完整分析产物

与 v3.12 的区别
----------------
1. **Lee 2012 进入合并**。v3.12 把它记为"无可提取效应量"，这是错的。其 Table 3 脚注原文：
   "Data are means ± SE. Values for the control group are imputed means ± SE. Change values
    for the intervention groups are the difference as compared with control with SE."
   即 105.7 / 37.4 是**相对对照的组间差值（net）**，连同 SE 可直接作效应量。
2. **按基线血糖状态分层**。Slentz 2009（原文排除糖尿病、基线 FPG 91–94 mg/dl、自述
   normal glucose tertiles）与 Lee 2012（nondiabetic，FPG 94.8–96.9）属正常血糖肥胖；
   DOSE-EX、Karstoft 属 T2D。v3.12 以"人群正常血糖"排除 HEARTY 却保留 Slentz，
   同模板双标；v3.13 改为分层全报，不再以血糖状态作排除理由。
3. **臂数 6 → 8，研究数 3 → 4**（DOSE-EX 2 臂、Karstoft 1、Slentz 3、Lee 2）。

Lee 的分泌轴取 first-phase insulin（急性相），与 Slentz 的 AIRg 同为动态急性指标；
second-phase 亦报告但未采用（敏感性见 analysis/v313_lee_secondphase.json）。

输出: analysis/meta_A_v313.json, analysis/decomposition_v313.csv, analysis/v313_summary.md

用法: python3 scripts/m2_v313.py
"""
import csv
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m2_meta as M

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Lee 2012: Table 1 基线 mean±SD；Table 3 ITT Δ vs control mean±SE
LEE = {
    "DI": {"base": {"AE": 538.9, "RE": 545.5}, "delta": {"AE": (105.7, 78.0), "RE": (37.4, 77.1)}},
    "SI": {"base": {"AE": 2.2, "RE": 2.9}, "delta": {"AE": (0.4, 0.2), "RE": (0.8, 0.2)}},
    "SEC": {"base": {"AE": 292.7, "RE": 213.9}, "delta": {"AE": (3.9, 22.4), "RE": (-38.2, 21.5)}},
}
LEE_SEC_SECOND = {"AE": (-36.2, 23.2, 344.0), "RE": (-48.6, 22.1, 241.4)}  # Δ, SE, baseline

LABEL = {  # 血糖分层标签
    "DOSE-EX 2023": "T2D", "Karstoft 2014": "T2D",
    "Slentz 2009": "normoglycaemic obesity", "Lee 2012": "normoglycaemic obesity",
}
ORDER = ["DOSE-EX 2023", "Karstoft 2014", "Slentz 2009", "Lee 2012"]


def lee_records():
    out = []
    for q, spec in LEE.items():
        for arm in ["AE", "RE"]:
            b = spec["base"][arm]
            d, se_d = spec["delta"][arm]
            post = b + d
            out.append({"study": "Lee 2012", "arm": arm, "q": q,
                        "ln": math.log(post / b), "se": se_d / post,
                        "src": "Table 3 (Δ vs control, ITT)"})
    return out


def build():
    rec = M.build_set_A() + lee_records()
    for r in rec:
        r["stratum"] = LABEL[r["study"]]
    rec.sort(key=lambda r: (ORDER.index(r["study"]), r["arm"], r["q"]))
    return rec


def pool(rows):
    return M.dl_kh([r["ln"] for r in rows], [r["se"] for r in rows])


def weights(rows):
    g = pool(rows)
    w = np.array([1 / (r["se"] ** 2 + g["tau2"]) for r in rows])
    return g, list(w / w.sum() * 100)


def main():
    rec = build()
    out = {"per_arm_data": [], "A_exercise_vs_control": {}, "strata": {}, "sensitivity": {}}

    # 权重按轴计算
    for q in ["DI", "SI", "SEC"]:
        rows = [r for r in rec if r["q"] == q]
        g, w = weights(rows)
        wmap = {(r["study"], r["arm"]): wi for r, wi in zip(rows, w)}
        for r in rows:
            out["per_arm_data"].append({
                "study": r["study"], "arm": r["arm"], "q": q, "ln": r["ln"], "se": r["se"],
                "pct": round((math.exp(r["ln"]) - 1) * 100, 1),
                "ci_lo_pct": round((math.exp(r["ln"] - 1.959964 * r["se"]) - 1) * 100, 1),
                "ci_hi_pct": round((math.exp(r["ln"] + 1.959964 * r["se"]) - 1) * 100, 1),
                "weight_pct": round(wmap[(r["study"], r["arm"])], 1),
                "stratum": r["stratum"], "src": r["src"],
            })
        out["A_exercise_vs_control"][q] = {
            "k": g["k"], "mu_ln": g["mu_ln"], "pct": round(g["pct"], 1),
            "ci_lo_pct": round(g["ci_lo_pct"], 1), "ci_hi_pct": round(g["ci_hi_pct"], 1),
            "p": float(f"{g['p']:.4g}"), "I2": round(g["I2"], 1), "tau2": round(g["tau2"], 5),
            "method": "DL tau2 + Knapp-Hartung",
            "arms": [r["study"] + "/" + r["arm"] for r in rows],
        }

    # 分层
    for sname, keep in [("T2D", {"DOSE-EX 2023", "Karstoft 2014"}),
                        ("normoglycaemic obesity", {"Slentz 2009", "Lee 2012"})]:
        out["strata"][sname] = {}
        for q in ["DI", "SI", "SEC"]:
            rows = [r for r in rec if r["q"] == q and r["study"] in keep]
            if len(rows) < 2:
                out["strata"][sname][q] = {"k": len(rows), "pooled": False}
                continue
            g = pool(rows)
            out["strata"][sname][q] = {
                "k": g["k"], "pooled": True, "pct": round(g["pct"], 1),
                "ci_lo_pct": round(g["ci_lo_pct"], 1), "ci_hi_pct": round(g["ci_hi_pct"], 1),
                "p": float(f"{g['p']:.4g}"), "I2": round(g["I2"], 1), "tau2": round(g["tau2"], 5),
                "arms": [r["study"] + "/" + r["arm"] for r in rows],
            }

    # 敏感性：每研究一臂（取该研究该轴权重最大臂）
    for q in ["DI", "SI", "SEC"]:
        rows = [r for r in rec if r["q"] == q]
        g, w = weights(rows)
        best = {}
        for r, wi in zip(rows, w):
            if r["study"] not in best or wi > best[r["study"]][1]:
                best[r["study"]] = (r, wi)
        picked = [v[0] for v in best.values()]
        if len(picked) >= 2:
            gs = pool(picked)
            out["sensitivity"][q] = {
                "k": gs["k"], "pct": round(gs["pct"], 1),
                "ci_lo_pct": round(gs["ci_lo_pct"], 1), "ci_hi_pct": round(gs["ci_hi_pct"], 1),
                "p": float(f"{gs['p']:.4g}"), "I2": round(gs["I2"], 1),
                "arms": [r["study"] + "/" + r["arm"] for r in picked],
            }

    os.makedirs(os.path.join(ROOT, "analysis"), exist_ok=True)
    with open(os.path.join(ROOT, "analysis/meta_A_v313.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    # 分解表（需三轴齐全的臂）
    rows = []
    by_arm = {}
    for r in rec:
        by_arm.setdefault((r["study"], r["arm"]), {})[r["q"]] = r
    for (study, arm), d in sorted(by_arm.items(), key=lambda kv: (ORDER.index(kv[0][0]), kv[0][1])):
        if len(d) == 3:
            rows.append({
                "study": study, "arm": arm, "stratum": LABEL[study],
                "dln_DI": round(d["DI"]["ln"], 4), "dln_SI": round(d["SI"]["ln"], 4),
                "dln_SEC": round(d["SEC"]["ln"], 4),
                "residual": round(d["DI"]["ln"] - d["SI"]["ln"] - d["SEC"]["ln"], 4),
                "si_share": round(d["SI"]["ln"] / d["DI"]["ln"], 2) if d["DI"]["ln"] else None,
                "sec_share": round(d["SEC"]["ln"] / d["DI"]["ln"], 2) if d["DI"]["ln"] else None,
            })
    with open(os.path.join(ROOT, "analysis/decomposition_v313.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # Lee second-phase 敏感性
    sp = {}
    for arm, (d, se, b) in LEE_SEC_SECOND.items():
        post = b + d
        sp[arm] = {"ln": round(math.log(post / b), 4), "se": round(se / post, 4)}
    with open(os.path.join(ROOT, "analysis/v313_lee_secondphase.json"), "w", encoding="utf-8") as f:
        json.dump(sp, f, indent=2)

    # 摘要打印
    L = ["# v3.13 analysis summary", ""]
    L.append("## Pooled, all arms")
    L.append("| Axis | k arms | k studies | Pooled | 95% CI | p | I² | τ² |")
    L.append("|---|---|---|---|---|---|---|---|")
    for q, name in [("DI", "Disposition index"), ("SI", "Insulin sensitivity"), ("SEC", "Insulin secretion")]:
        m = out["A_exercise_vs_control"][q]
        ns = len({a.split("/")[0] for a in m["arms"]})
        L.append(f"| {name} | {m['k']} | {ns} | {m['pct']:+.1f}% | {m['ci_lo_pct']:+.1f} to {m['ci_hi_pct']:+.1f} | "
                 f"{m['p']} | {m['I2']}% | {m['tau2']} |")
    L.append("")
    L.append("## By stratum")
    L.append("| Stratum | Axis | k | Pooled | 95% CI | p | I² |")
    L.append("|---|---|---|---|---|---|---|")
    for s in ["T2D", "normoglycaemic obesity"]:
        for q in ["DI", "SI", "SEC"]:
            d = out["strata"][s][q]
            if not d["pooled"]:
                L.append(f"| {s} | {q} | {d['k']} | not pooled | | | |")
            else:
                L.append(f"| {s} | {q} | {d['k']} | {d['pct']:+.1f}% | {d['ci_lo_pct']:+.1f} to {d['ci_hi_pct']:+.1f} | {d['p']} | {d['I2']}% |")
    L.append("")
    L.append("## One arm per study")
    for q in ["DI", "SI", "SEC"]:
        if q in out["sensitivity"]:
            d = out["sensitivity"][q]
            L.append(f"- {q}: {d['pct']:+.1f}% ({d['ci_lo_pct']:+.1f}, {d['ci_hi_pct']:+.1f}); p={d['p']}; I²={d['I2']}%; k={d['k']} — {', '.join(d['arms'])}")
    L.append("")
    L.append("## Decomposition (arms reporting all three axes)")
    L.append("| Study/arm | Δln DI | Δln Si | Δln Sec | Residual |")
    L.append("|---|---|---|---|---|")
    for r in rows:
        L.append(f"| {r['study']} {r['arm']} | {r['dln_DI']:+.4f} | {r['dln_SI']:+.4f} | {r['dln_SEC']:+.4f} | {r['residual']:+.4f} |")
    with open(os.path.join(ROOT, "analysis/v313_summary.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
