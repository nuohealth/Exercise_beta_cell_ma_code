#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v5.0 图件的单一数据源。

所有图（Fig. 2 / Fig. 3）都从这里取逐臂数据，杜绝"正文 k=9、图还是 8 臂"的漂移。
数值口径与 scripts/a1p_add_study.py 完全一致（同一批 ln / SE），并逐项自检。

输出: analysis/arms_v500.json
"""
import json
import math
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ---- 原始 8 臂（v4.0 已发表口径，一个字节不改） ----
ORIG = {
    "DI": {
        ("DOSE-EX 2023", "HED"): (0.8628899551470398, 0.16085335950863347),
        ("DOSE-EX 2023", "MED"): (0.7178397931503168, 0.16274808236086336),
        ("Karstoft 2014", "IWT"): (0.5080216964332565, 0.13116726835138387),
        ("Slentz 2009", "LOWMOD"): (0.5026162977469957, 0.13552476903859137),
        ("Slentz 2009", "LOWVIG"): (0.19637502284164818, 0.108532266286991),
        ("Slentz 2009", "HIGHVIG"): (0.22098672646679193, 0.0957828570939168),
    },
    "SI": {
        ("DOSE-EX 2023", "HED"): (0.6043159668533296, 0.1561707556940531),
        ("DOSE-EX 2023", "MED"): (0.4054651081081644, 0.1575640236381383),
        ("Karstoft 2014", "IWT"): (0.40413088509502776, 0.09746328437917222),
        ("Slentz 2009", "LOWMOD"): (0.542574322080571, 0.09889769019129266),
        ("Slentz 2009", "LOWVIG"): (0.3011050927839216, 0.08820599358380829),
        ("Slentz 2009", "HIGHVIG"): (0.40082472855166235, 0.09353263171956532),
    },
    "SEC": {
        ("DOSE-EX 2023", "HED"): (0.2468600779315258, 0.06360982235087835),
        ("DOSE-EX 2023", "MED"): (0.3220834991691132, 0.06434423301016023),
        ("Slentz 2009", "LOWMOD"): (-0.06597680233078876, 0.08555408044437353),
        ("Slentz 2009", "LOWVIG"): (-0.14536883215782767, 0.09234497349850307),
        ("Slentz 2009", "HIGHVIG"): (-0.20908844189875272, 0.08720418418205479),
    },
}

# ---- Lee 2012 两臂：由基线与变化量按 §2.6 口径导出 ----
LEE = {
    "DI": {"base": {"AE": 538.9, "RE": 545.5},
           "delta": {"AE": (105.7, 78.0), "RE": (37.4, 77.1)}},
    "SI": {"base": {"AE": 2.2, "RE": 2.9},
           "delta": {"AE": (0.4, 0.2), "RE": (0.8, 0.2)}},
    "SEC": {"base": {"AE": 292.7, "RE": 213.9},
            "delta": {"AE": (3.9, 22.4), "RE": (-38.2, 21.5)}},
}

STRATUM = {
    "DOSE-EX 2023": "T2D",
    "Karstoft 2014": "T2D",
    "Slentz 2009": "NGT",
    "Lee 2012": "NGT",
    "Carr 2005": "IGT",
}

# ---- Carr 2005 (15677490)：Si 4.7±0.5 vs 3.3±0.3；AIRg 294.9±50.0 vs 305.4±30.0 ----
m1, s1 = 4.7, 0.5
m2, s2 = 3.3, 0.3
ln_Si = math.log(m1 / m2)
se_Si = math.sqrt((s1 / m1) ** 2 + (s2 / m2) ** 2)

a1m, a1s = 294.9, 50.0
a2m, a2s = 305.4, 30.0
ln_Ag = math.log(a1m / a2m)
se_Ag = math.sqrt((a1s / a1m) ** 2 + (a2s / a2m) ** 2)

CARR = {
    "SI": (ln_Si, se_Si),
    "SEC": (ln_Ag, se_Ag),
    "DI": (ln_Si + ln_Ag, math.sqrt(se_Si ** 2 + se_Ag ** 2)),
}


def lee_arm(q, arm):
    b = LEE[q]["base"][arm]
    d, se_d = LEE[q]["delta"][arm]
    post = b + d
    return (math.log(post / b), se_d / post)


def build():
    out = {}
    for q in ("DI", "SI", "SEC"):
        rows = []
        for (study, arm), (y, se) in ORIG[q].items():
            rows.append(dict(study=study, arm=arm, q=q, ln=y, se=se,
                             stratum=STRATUM[study], source="measured_original"))
        for arm in ("AE", "RE"):
            y, se = lee_arm(q, arm)
            rows.append(dict(study="Lee 2012", arm=arm, q=q, ln=y, se=se,
                             stratum=STRATUM["Lee 2012"], source="derived_baseline_delta"))
        y, se = CARR[q]
        rows.append(dict(study="Carr 2005", arm="", q=q, ln=y, se=se,
                         stratum=STRATUM["Carr 2005"],
                         source="derived_product" if q == "DI" else "derived_ratio"))
        out[q] = rows
    return out


def dl_kh(y, v):
    k = len(y)
    w = [1.0 / x for x in v]
    sw = sum(w)
    sw2 = sum(x * x for x in w)
    mu_fe = sum(wi * yi for wi, yi in zip(w, y)) / sw
    Q = sum(wi * (yi - mu_fe) ** 2 for wi, yi in zip(w, y))
    C = sw - sw2 / sw
    tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    wr = [1.0 / (vi + tau2) for vi in v]
    swr = sum(wr)
    mu = sum(wi * yi for wi, yi in zip(wr, y)) / swr
    Qre = sum(wi * (yi - mu) ** 2 for wi, yi in zip(wr, y))
    var_hk = Qre / ((k - 1) * swr) if k > 1 else 1.0 / swr
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    return dict(k=k, mu=mu, se=math.sqrt(var_hk), tau2=tau2, I2=I2, Q=Q, df=k - 1,
                pct=(math.exp(mu) - 1) * 100)


if __name__ == "__main__":
    arms = build()
    # 自检：与 meta_A_a1p.json 完全一致，否则拒绝出图
    ref = json.load(open(os.path.join(BASE, "analysis/meta_A_a1p.json"), encoding="utf-8"))
    ok = True
    for q in ("DI", "SI", "SEC"):
        rows = arms[q]
        d = dl_kh([r["ln"] for r in rows], [r["se"] ** 2 for r in rows])
        r = ref[q]["new"]
        # Knapp-Hartung interval uses the t critical value at df = k-1, not a
        # fixed 1.96: at k=9 that is 2.306, at k=8 it is 2.365.
        from scipy.stats import t as student_t
        tcv = student_t.ppf(0.975, d["df"])
        lo = round(100 * (math.exp(d["mu"] - tcv * d["se"]) - 1), 1)
        hi = round(100 * (math.exp(d["mu"] + tcv * d["se"]) - 1), 1)
        same = (d["k"] == r["k"]
                and abs(round(d["pct"], 1) - r["pct"]) < 0.05
                and abs(lo - r["ci"][0]) < 0.25
                and abs(hi - r["ci"][1]) < 0.25)
        print(f"{q}: k={d['k']} pct={round(d['pct'],1)} "
              f"(ref k={r['k']} pct={r['pct']})  {'OK' if same else 'MISMATCH'}")
        ok = ok and same
    if not ok:
        raise SystemExit("图件数据源与正文不一致，拒绝出图")
    path = os.path.join(BASE, "analysis", "arms_v500.json")
    json.dump(arms, open(path, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print("→", path)
