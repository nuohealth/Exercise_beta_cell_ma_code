#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
m4_interaction_test.py — 分层 × 分泌 交互作用检验（审稿人前置防御）
=====================================================================
背景：v3.14 §3.3/§4.5 主张"两层分泌方向相反（T2D +32.8% vs 正常血糖 −10.0%）"，
      但未报告交互作用检验。第二轮外审指出"opposition is significant"缺统计支持。
      本脚本补算该检验，判断该主张是否站得住。

单一数据源：analysis/meta_A_v313.json（per_arm_data 中 q == "SEC"）
输出：analysis/interaction_test_v314.json

三种互补检验（均在同一数据源上）：
  M1  混合效应亚组差值检验（Borenstein 亚组分析）：diff = mu_T2D − mu_NGT，
      se = sqrt(se_T2D² + se_NGT²)；同时给正态 z 与 Knapp–Hartung t（df = k1+k2−2）
  M2  固定效应 Q-between：Q_b = Q_total − Q_within，df = 层数−1
  M3  稳健性：剔除 Lee 2012（分层标签为推断而非 OGTT/HbA1c 实测）后重算 NGT 层

注意：T2D 分泌层仅 k=2（DOSE-EX MED/HED，共享同一 CON 对照），
      故 M1 的 KH t 版本 df=5 偏乐观；M2 与 M3 用于交叉印证方向一致性。

口径与 m2_v313.py 完全一致：DL τ² + Knapp–Hartung 调整。
"""
import json
import math
import os

import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "analysis", "meta_A_v313.json")
OUT = os.path.join(ROOT, "analysis", "interaction_test_v314.json")


def dl_kh(y, se):
    """DerSimonian–Laird τ² + Knapp–Hartung 调整（与 m2_meta.dl_kh 同口径）"""
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
    q_kh = (wr * (y - mu) ** 2).sum() / df if df > 0 else 1.0
    return mu, se_mu * math.sqrt(max(q_kh, 1e-9)), tau2, float(Q), df


def main():
    data = json.load(open(SRC, encoding="utf-8"))
    rows = [r for r in data["per_arm_data"] if r["q"] == "SEC"]

    # 分层合并
    grp = {}
    for r in rows:
        grp.setdefault(r["stratum"], []).append((r["ln"], r["se"]))

    strata = {}
    for name, v in grp.items():
        mu, se, tau2, Q, df = dl_kh([a for a, _ in v], [b for _, b in v])
        strata[name] = {
            "k": len(v), "mu_ln": round(mu, 6), "se": round(se, 6),
            "pct": round((math.exp(mu) - 1) * 100, 1),
            "tau2": round(tau2, 6), "Q": round(Q, 3), "df": df,
            "arms": [f"{r['study']}/{r['arm']}" for r in rows if r["stratum"] == name],
        }
        print(f"{name:<26s} k={len(v)}  {100*(math.exp(mu)-1):+6.1f}%  se={se:.4f}  tau2={tau2:.5f}")

    a, sa, ka = strata["T2D"]["mu_ln"], strata["T2D"]["se"], strata["T2D"]["k"]
    b, sb, kb = (strata["normoglycaemic obesity"]["mu_ln"],
                 strata["normoglycaemic obesity"]["se"],
                 strata["normoglycaemic obesity"]["k"])

    # M1 亚组差值
    diff = a - b
    se_diff = math.sqrt(sa**2 + sb**2)
    z = diff / se_diff
    df_t = ka + kb - 2
    m1 = {
        "diff_ln": round(diff, 6),
        "ratio_pct": round((math.exp(diff) - 1) * 100, 1),
        "se_diff": round(se_diff, 6),
        "z": round(z, 3),
        "p_normal": round(float(2 * (1 - stats.norm.cdf(abs(z)))), 6),
        "p_kh_t": round(float(2 * (1 - stats.t.cdf(abs(z), df_t))), 6),
        "df_t": df_t,
        "note": "T2D 层 k=2 共享同一对照，KH t 的 df 偏乐观；以 M2/M3 交叉印证",
    }
    print(f"\n[M1 亚组差值] diff_ln={diff:+.4f} ({m1['ratio_pct']:+.1f}%)  "
          f"z={z:.3f}  p(normal)={m1['p_normal']:.2e}  p(KH t, df={df_t})={m1['p_kh_t']:.4f}")

    # M2 Q-between（固定效应）
    y = np.array([r["ln"] for r in rows])
    se = np.array([r["se"] for r in rows])
    st = [r["stratum"] for r in rows]
    w = 1 / se**2
    fe = (w * y).sum() / w.sum()
    Qtot = float((w * (y - fe) ** 2).sum())
    Qw = 0.0
    for s in set(st):
        m = np.array([t == s for t in st])
        ws, ys = w[m], y[m]
        mu = (ws * ys).sum() / ws.sum()
        Qw += float((ws * (ys - mu) ** 2).sum())
    Qb = Qtot - Qw
    dfb = len(set(st)) - 1
    m2 = {
        "Q_total": round(Qtot, 3), "Q_within": round(Qw, 3),
        "Q_between": round(Qb, 3), "df": dfb,
        "p": round(float(1 - stats.chi2.cdf(Qb, dfb)), 8),
    }
    print(f"[M2 Q-between] Q_b={Qb:.3f}  df={dfb}  p={m2['p']:.2e}")

    # M3 剔除 Lee 2012 后的 NGT 层
    ngt_wo = [(r["ln"], r["se"]) for r in rows
              if r["stratum"] == "normoglycaemic obesity" and r["study"] != "Lee 2012"]
    mu3, se3, tau3, Q3, df3 = dl_kh([x for x, _ in ngt_wo], [x for _, x in ngt_wo])
    d3 = a - mu3
    s3 = math.sqrt(sa**2 + se3**2)
    m3 = {
        "k": len(ngt_wo), "pct": round((math.exp(mu3) - 1) * 100, 1),
        "mu_ln": round(mu3, 6), "se": round(se3, 6), "tau2": round(tau3, 6),
        "diff_vs_T2D_ln": round(d3, 6),
        "p_normal": round(float(2 * (1 - stats.norm.cdf(abs(d3 / s3)))), 6),
        "arms": [f"{r['study']}/{r['arm']}" for r in rows
                 if r["stratum"] == "normoglycaemic obesity" and r["study"] != "Lee 2012"],
    }
    print(f"[M3 去 Lee 后 NGT] k={len(ngt_wo)}  {m3['pct']:+.1f}%  "
          f"diff vs T2D={d3:+.4f}  p={m3['p_normal']:.2e}")

    verdict = (m1["p_normal"] < 0.05 and m2["p"] < 0.05 and m3["p_normal"] < 0.05)
    out = {
        "question": "两层分泌效应方向相反，其差异是否统计显著（交互作用检验）",
        "strata": strata,
        "M1_subgroup_difference": m1,
        "M2_Q_between": m2,
        "M3_sensitivity_excluding_Lee": m3,
        "verdict_opposition_significant": verdict,
        "interpretation": (
            "三种检验一致显著，故‘opposition’有统计支持；v3.14 未报告该检验是呈现缺口，"
            "而非结论缺陷。建议在 §3.3 补一句交互检验结果，并在 §4.5 将"
            "‘what is not compatible with no effect is their opposition’改为引用该检验的表述。"
            if verdict else
            "检验未达显著，应弱化‘opposition’表述为‘两层点估计方向相反但均不显著’。"
        ),
        "source": "analysis/meta_A_v313.json :: per_arm_data[q=='SEC']",
    }
    json.dump(out, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print(f"\n判定 opposition 显著 = {verdict}")
    print(f"已写出 {OUT}")


if __name__ == "__main__":
    main()
