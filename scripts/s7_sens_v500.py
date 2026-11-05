#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S7 v5.0：k=9 全套敏感性分析。

单一数据源: analysis/arms_v500.json（自身已对 meta_A_a1p.json 自校验）
估计量: DerSimonian–Laird τ² + Knapp–Hartung 区间，与正文 §2.6 同一实现。

输出: analysis/s7_sensitivity_v500.json
"""
import json
import math
import os

from scipy.stats import t as student_t

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AXES = ("DI", "SI", "SEC")


def dl_kh(y, v):
    k = len(y)
    w = [1.0 / x for x in v]
    sw = sum(w)
    mu_fe = sum(wi * yi for wi, yi in zip(w, y)) / sw
    Q = sum(wi * (yi - mu_fe) ** 2 for wi, yi in zip(w, y))
    C = sw - sum(x * x for x in w) / sw
    tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 and k > 1 else 0.0
    wr = [1.0 / (vi + tau2) for vi in v]
    swr = sum(wr)
    mu = sum(wi * yi for wi, yi in zip(wr, y)) / swr
    Qre = sum(wi * (yi - mu) ** 2 for wi, yi in zip(wr, y))
    df = k - 1
    var = Qre / (df * swr) if df > 0 else 1.0 / swr
    se = math.sqrt(var)
    tcv = student_t.ppf(0.975, df) if df > 0 else 0.0
    p = 2 * (1 - student_t.cdf(abs(mu / se), df)) if df > 0 else float("nan")
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 and df > 0 else 0.0
    return dict(k=k, pct=round((math.exp(mu) - 1) * 100, 1),
                ci=[round((math.exp(mu - tcv * se) - 1) * 100, 1),
                    round((math.exp(mu + tcv * se) - 1) * 100, 1)],
                p=round(p, 4), I2=round(I2, 1), tau2=round(tau2, 3),
                mu=mu, se=se, df=df, Q=round(Q, 2))


def heaviest(rows):
    """每研究取权重最大的一臂（权重 = 1/(se^2+tau2)，用 τ²=0 的近似排序，与旧版一致）。"""
    best = {}
    for r in rows:
        w = 1.0 / r["se"] ** 2
        if r["study"] not in best or w > best[r["study"]][0]:
            best[r["study"]] = (w, r)
    return [v[1] for v in best.values()]


def main():
    arms = json.load(open(os.path.join(BASE, "analysis/arms_v500.json"), encoding="utf-8"))
    ref = json.load(open(os.path.join(BASE, "analysis/meta_A_a1p.json"), encoding="utf-8"))
    out = {}

    for q in AXES:
        rows = arms[q]
        res = {}
        res["main"] = dl_kh([r["ln"] for r in rows], [r["se"] ** 2 for r in rows])
        # 闸门：主分析必须与正文一致
        m, r0 = res["main"], ref[q]["new"]
        assert m["k"] == r0["k"] and abs(m["pct"] - r0["pct"]) < 0.05, (q, m, r0)

        op = heaviest(rows)
        res["one_per_study"] = dl_kh([r["ln"] for r in op], [r["se"] ** 2 for r in op])
        res["one_per_study_arms"] = [f"{r['study']}/{r['arm'] or '-'}" for r in op]

        for study, tag in [("DOSE-EX 2023", "excl_DOSEEX"), ("Lee 2012", "excl_Lee"),
                           ("Slentz 2009", "excl_Slentz"), ("Karstoft 2014", "excl_Karstoft"),
                           ("Carr 2005", "excl_Carr")]:
            sub = [r for r in rows if r["study"] != study]
            res[tag] = dl_kh([r["ln"] for r in sub], [r["se"] ** 2 for r in sub])
        out[q] = res

    # 层级权重份额（DI 轴，用于 D24/D26 的 65.5% 口径 → 现需按三层重算）
    di = arms["DI"]
    shares = {}
    for st in ("T2D", "NGT", "IGT"):
        sub = [r for r in di if r["stratum"] == st]
        m = dl_kh([r["ln"] for r in sub], [r["se"] ** 2 for r in sub])
        tau2 = m["tau2"]
        wr = [1.0 / (r["se"] ** 2 + tau2) for r in sub]
        swr = sum(1.0 / (r["se"] ** 2 + tau2) for r in di)
        shares[st] = dict(arms=len(sub), studies=len({r["study"] for r in sub}),
                          weight=round(100 * sum(wr) / swr, 1))

    path = os.path.join(BASE, "analysis", "s7_sensitivity_v500.json")
    json.dump(dict(axes=out, di_stratum_weight=shares),
              open(path, "w", encoding="utf-8"), indent=1, ensure_ascii=False)

    for q in AXES:
        print(f"== {q}")
        for k, v in out[q].items():
            if k == "one_per_study_arms":
                print("   arms:", ", ".join(v)); continue
            print(f"   {k:16s} k={v['k']:2d} {v['pct']:+7.1f}% "
                  f"({v['ci'][0]:+7.1f},{v['ci'][1]:+7.1f}) p={v['p']:<8} I2={v['I2']}")
    print("\nDI 轴层级权重:", shares)
    print("→", path)


if __name__ == "__main__":
    main()
