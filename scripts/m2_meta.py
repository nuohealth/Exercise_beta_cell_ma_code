#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
随机效应 Meta 分析：ln(DI)、ln(Si)、ln(Sec) 三层分别合并
方法：DerSimonian-Laird τ² + Knapp-Hartung 校正 CI（k 小，必须用 KH）
两套对比：
  A. 运动 vs 无干预对照        —— 回答"运动有没有用"
  B. 等量减重前提下加运动      —— 回答"运动的作用是否独立于减重"（核心）
输出：analysis/meta_A.json, analysis/meta_B.json, analysis/meta_table.csv
"""
import json, math, os, csv
import numpy as np
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "analysis")
os.makedirs(OUT, exist_ok=True)
Z = 1.959964


def ln_se(lo, hi):
    return (math.log(hi) - math.log(lo)) / (2 * Z)


def dl_kh(y, se):
    """DerSimonian-Laird τ² + Knapp-Hartung 调整"""
    y, se = np.asarray(y, float), np.asarray(se, float)
    k = len(y)
    w = 1 / se**2
    fe = (w * y).sum() / w.sum()
    Q = (w * (y - fe) ** 2).sum()
    df = k - 1
    C = w.sum() - (w**2).sum() / w.sum()
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    wr = 1 / (se**2 + tau2)
    mu = (wr * y).sum() / wr.sum()
    se_mu = math.sqrt(1 / wr.sum())
    # Knapp-Hartung
    q_kh = (wr * (y - mu) ** 2).sum() / df if df > 0 else 1.0
    se_kh = se_mu * math.sqrt(max(q_kh, 1e-9))
    tcrit = stats.t.ppf(0.975, df) if df > 0 else Z
    p = 2 * (1 - stats.t.cdf(abs(mu / se_kh), df)) if df > 0 else 2 * (1 - stats.norm.cdf(abs(mu / se_kh)))
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    return {
        "k": k, "mu_ln": round(mu, 4), "se": round(se_kh, 4),
        "pct": round((math.exp(mu) - 1) * 100, 1),
        "ci_lo_pct": round((math.exp(mu - tcrit * se_kh) - 1) * 100, 1),
        "ci_hi_pct": round((math.exp(mu + tcrit * se_kh) - 1) * 100, 1),
        "p": round(float(p), 4), "tau2": round(tau2, 4), "I2": round(I2, 1),
        "Q": round(float(Q), 3), "df": df,
        "method": "DL τ² + Knapp-Hartung",
    }


# ══════════════════════════════════════════════════════════════════
# 数据：全部来自原文表格，逐条注明出处
# ══════════════════════════════════════════════════════════════════

# DOSE-EX 37127822 Nature Metab 2023, Table 2（相对 CON 的几何均值比 %，95%CI）
DOSEEX_VS_CON = {
    "HED": {"DI": (137, 73, 225), "SI": (83, 35, 149), "SEC": (28, 13, 45)},
    "MED": {"DI": (105, 49, 182), "SI": (50, 10, 104), "SEC": (38, 22, 57)},
    "DCON": {"DI": (58, 16, 116), "SI": (18, -13, 60), "SEC": (33, 18, 50)},
}
# 等量减重上加运动（vs DCON）
DOSEEX_VS_DCON = {
    "HED": {"DI": (50, 10, 104), "SI": (55, 15, 109), "SEC": (-4, -15, 9)},
    "MED": {"DI": (29, -5, 77), "SI": (27, -7, 71), "SEC": (4, -8, 18)},
}

# Karstoft 25099941 Diabetologia 2014（IWT vs CON，%±SEM）
KARSTOFT = {
    "DI": (66.2, 21.8),
    "SI": (49.8, 14.6),
    "SEC": (0.0, None),  # 原文：ISR 所有阶段无变化；p 值未给，SE 不可得
}

# Slentz 19592624 Diabetes Care 2009, Table 2（算术 mean±SD, baseline & change）
SLENTZ = {
    # arm: (DI b, DI chg, DI chgSD, AIRg b, AIRg chg, AIRg chgSD, Si b, Si chg, Si chgSD, n)
    "CON":     (1562, -65, 893, 490.6, 22, 233, 3.5, -0.35, 1.7, 58),
    "LOWMOD":  (1270, 742, 1680, 457.3, -10, 207, 3.1, 1.7, 2.5, 57),
    "LOWVIG":  (1533, 255, 1023, 445.5, -43, 216, 3.7, 0.8, 1.8, 58),
    "HIGHVIG": (1305, 255, 688, 505.6, -77, 218, 3.2, 1.1, 2.1, 64),
}


def slentz_ln(b, chg, sd, n):
    """算术 mean±SD → ln(post/base) 及 delta-method SE"""
    post = b + chg
    if b <= 0 or post <= 0:
        return None, None
    lnr = math.log(post / b)
    se = sd / (post * math.sqrt(n))
    return lnr, se


def build_set_A():
    """运动 vs 无干预对照"""
    rec = []
    for arm in ["HED", "MED"]:
        for q in ["DI", "SI", "SEC"]:
            p, lo, hi = DOSEEX_VS_CON[arm][q]
            rec.append({"study": "DOSE-EX 2023", "arm": arm, "q": q,
                        "ln": math.log(1 + p / 100), "se": ln_se(1 + lo / 100, 1 + hi / 100),
                        "src": "Table 2 vs CON"})
    # Karstoft
    for q in ["DI", "SI"]:
        p, sem = KARSTOFT[q]
        rec.append({"study": "Karstoft 2014", "arm": "IWT", "q": q,
                    "ln": math.log(1 + p / 100), "se": (sem / 100) / (1 + p / 100),
                    "src": "text %±SEM"})
    # Slentz：臂内 ln 比值 − CON 臂内 ln 比值
    con = {}
    for q, idx in [("DI", (0, 1, 2)), ("SI", (6, 7, 8)), ("SEC", (3, 4, 5))]:
        b, chg, sd = SLENTZ["CON"][idx[0]], SLENTZ["CON"][idx[1]], SLENTZ["CON"][idx[2]]
        con[q] = slentz_ln(b, chg, sd, SLENTZ["CON"][9])
    for arm in ["LOWMOD", "LOWVIG", "HIGHVIG"]:
        v = SLENTZ[arm]
        for q, idx in [("DI", (0, 1, 2)), ("SI", (6, 7, 8)), ("SEC", (3, 4, 5))]:
            b, chg, sd = v[idx[0]], v[idx[1]], v[idx[2]]
            lnr, se = slentz_ln(b, chg, sd, v[9])
            if lnr is None:
                continue
            cln, cse = con[q]
            net = lnr - cln
            netse = math.sqrt(se**2 + cse**2)
            rec.append({"study": "Slentz 2009", "arm": arm, "q": q,
                        "ln": net, "se": netse, "src": "Table 2 (arm − CON)"})
    return rec


def main():
    A = build_set_A()
    res_A = {}
    for q in ["DI", "SI", "SEC"]:
        sub = [r for r in A if r["q"] == q]
        res_A[q] = dl_kh([r["ln"] for r in sub], [r["se"] for r in sub])
        res_A[q]["arms"] = [f"{r['study']}/{r['arm']}" for r in sub]

    # 研究层敏感性：每研究只取一个臂（DOSE-EX 取 HED、Slentz 取 LOWMOD）
    sens = {}
    for q in ["DI", "SI", "SEC"]:
        keep = [r for r in A if r["q"] == q and not (
            (r["study"] == "DOSE-EX 2023" and r["arm"] == "MED") or
            (r["study"] == "Slentz 2009" and r["arm"] in ("LOWVIG", "HIGHVIG")))]
        if len(keep) >= 2:
            sens[q] = dl_kh([r["ln"] for r in keep], [r["se"] for r in keep])

    # B 套：等量减重下加运动（单研究，不合并，直接报原文）
    B = {"note": "仅 DOSE-EX 一项研究提供该对比，不做合并；原文 Table 2 直接读取",
         "HED_vs_DCON": {}, "MED_vs_DCON": {}}
    for arm in ["HED", "MED"]:
        for q in ["DI", "SI", "SEC"]:
            p, lo, hi = DOSEEX_VS_DCON[arm][q]
            B[f"{arm}_vs_DCON"][q] = {"pct": p, "ci": [lo, hi],
                                      "signif": not (lo <= 0 <= hi)}

    json.dump({"A_exercise_vs_control": res_A, "A_sensitivity_one_arm_per_study": sens,
               "per_arm_data": A},
              open(os.path.join(OUT, "meta_A.json"), "w"), indent=2, ensure_ascii=False)
    json.dump(B, open(os.path.join(OUT, "meta_B.json"), "w"), indent=2, ensure_ascii=False)

    # 打印
    print("=" * 84)
    print("A. 运动 vs 无干预对照（DL τ² + Knapp-Hartung）")
    print("=" * 84)
    print(f"{'量':<10}{'k':>3}{'合并%':>9}{'95%CI':>22}{'p':>9}{'I²':>8}{'τ²':>9}")
    for q, lab in [("DI", "DI"), ("SI", "敏感性Si"), ("SEC", "分泌Sec")]:
        r = res_A[q]
        ci = "[%.1f, %.1f]" % (r["ci_lo_pct"], r["ci_hi_pct"])
        print(f"{lab:<10}{r['k']:>3}{r['pct']:>9.1f}{ci:>22}"
              f"{r['p']:>9.3f}{r['I2']:>8.1f}{r['tau2']:>9.3f}")
    print("\n研究层敏感性（每研究 1 臂）:")
    for q in ["DI", "SI", "SEC"]:
        if q in sens:
            r = sens[q]
            print(f"  {q:<5}k={r['k']}  {r['pct']:+.1f}% [{r['ci_lo_pct']:+.1f}, {r['ci_hi_pct']:+.1f}] p={r['p']:.3f} I²={r['I2']:.1f}")
    print("\n" + "=" * 84)
    print("B. 等量减重前提下加运动（DOSE-EX vs DCON，单研究）")
    print("=" * 84)
    for arm in ["HED", "MED"]:
        line = f"  {arm} vs DCON: "
        for q in ["DI", "SI", "SEC"]:
            d = B[f"{arm}_vs_DCON"][q]
            line += f"{q} {d['pct']:+.0f}% [{d['ci'][0]:+.0f},{d['ci'][1]:+.0f}]{'*' if d['signif'] else ''}  "
        print(line)
    print("\n  * = 95%CI 不含 0")
    print("  → 等量减重下加运动：DI 显著上升、敏感性显著上升、**分泌 0 变化**（HED −4%，CI −15~+9）")


if __name__ == "__main__":
    main()
