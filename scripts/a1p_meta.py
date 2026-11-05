"""
A1' 执行：纳入 15677490，重跑 DI / Si / Sec 三轴
第一步：先复现 v4.0 的 meta_A.json（不通过则说明实现不对，禁止继续）
第二步：加入 15677490 后重算
"""
import json, math

# ---------- 工具函数 ----------
def dl_kh(y, v):
    """DL tau^2 + Knapp-Hartung，与原稿 §2.6 及 meta_A.json 口径一致"""
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
    # HK
    var_hk = Qre / ((k - 1) * swr) if k > 1 else (1.0 / swr)
    se = math.sqrt(var_hk)
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    return dict(k=k, mu_ln=mu, se=se, tau2=tau2, I2=I2, Q=Q, df=k - 1,
                pct=(math.exp(mu) - 1) * 100)


def t_crit(df, alpha=0.05):
    # 双尾 t 临界值（小 df 精确表，df<=30）
    T = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365,
         8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160,
         14: 2.145, 15: 2.131, 16: 2.120, 17: 2.110, 18: 2.101, 19: 2.093,
         20: 2.086, 21: 2.080, 22: 2.074, 23: 2.069, 24: 2.064, 25: 2.060,
         26: 2.056, 27: 2.052, 28: 2.048, 29: 2.045, 30: 2.042}
    return T.get(df, 1.96)


def p_from_t(t, df):
    # 双尾 p，用不完全 beta 的粗算（足够用于报告）
    x = df / (df + t * t)
    # 正则不完全 beta I_x(df/2, 1/2)
    a, b = df / 2.0, 0.5
    if x <= 0: return 1.0
    if x >= 1: return 0.0
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    def betacf(aa, bb, xx):
        MAXIT, EPS, FPMIN = 200, 3e-16, 1e-300
        qab, qap, qam = aa + bb, aa + 1, aa - 1
        c = 1.0; d = 1 - qab * xx / qap
        if abs(d) < FPMIN: d = FPMIN
        d = 1 / d; h = d
        for m in range(1, MAXIT + 1):
            m2 = 2 * m
            aa2 = m * (bb - m) * xx / ((qam + m2) * (aa + m2))
            d = 1 + aa2 * d
            if abs(d) < FPMIN: d = FPMIN
            c = 1 + aa2 / c
            if abs(c) < FPMIN: c = FPMIN
            d = 1 / d; h *= d * c
            aa2 = -(aa + m) * (qab + m) * xx / ((aa + m2) * (qap + m2))
            d = 1 + aa2 * d
            if abs(d) < FPMIN: d = FPMIN
            c = 1 + aa2 / c
            if abs(c) < FPMIN: c = FPMIN
            d = 1 / d; de = d * c; h *= de
            if abs(de - 1) < EPS: break
        return h
    bt = math.exp(-lbeta + a * math.log(x) + b * math.log(1 - x))
    if x < (a + 1) / (a + b + 2):
        return bt * betacf(a, b, x) / a
    return 1 - bt * betacf(b, a, 1 - x) / b


def summarize(d, label):
    df = d["df"]
    tc = t_crit(df)
    lo = math.exp(d["mu_ln"] - tc * d["se"])
    hi = math.exp(d["mu_ln"] + tc * d["se"])
    tval = d["mu_ln"] / d["se"]
    p = p_from_t(abs(tval), df)
    return dict(label=label, k=d["k"], pct=round(d["pct"], 1),
                ci=(round((lo - 1) * 100, 1), round((hi - 1) * 100, 1)),
                p=round(p, 4), tau2=round(d["tau2"], 4), I2=round(d["I2"], 1),
                Q=round(d["Q"], 3), mu=round(d["mu_ln"], 4), se=round(d["se"], 4))


# ---------- 第一步：复现 v4.0 ----------
BASE = {
    "DI": [
        ("DOSE-EX/HED", 0.8628899551470398, 0.16085335950863347),
        ("DOSE-EX/MED", 0.7178397931503168, 0.16274808236086336),
        ("Karstoft/IWT", 0.5080216964332565, 0.13116726835138387),
        ("Slentz/LOWMOD", 0.5026162977469957, 0.13552476903859137),
        ("Slentz/LOWVIG", 0.19637502284164818, 0.108532266286991),
        ("Slentz/HIGHVIG", 0.22098672646679193, 0.0957828570939168),
    ],
    "SI": [
        ("DOSE-EX/HED", 0.6043159668533296, 0.1561707556940531),
        ("DOSE-EX/MED", 0.4054651081081644, 0.1575640236381383),
        ("Karstoft/IWT", 0.40413088509502776, 0.09746328437917222),
        ("Slentz/LOWMOD", 0.542574322080571, 0.09889769019129266),
        ("Slentz/LOWVIG", 0.3011050927839216, 0.08820599358380829),
        ("Slentz/HIGHVIG", 0.40082472855166235, 0.09353263171956532),
    ],
    "SEC": [
        ("DOSE-EX/HED", 0.2468600779315258, 0.06360982235087835),
        ("DOSE-EX/MED", 0.3220834991691132, 0.06434423301016023),
        ("Slentz/LOWMOD", -0.06597680233078876, 0.08555408044437353),
        ("Slentz/LOWVIG", -0.14536883215782767, 0.09234497349850307),
        ("Slentz/HIGHVIG", -0.20908844189875272, 0.08720418418205479),
    ],
}

EXPECT = {
    "DI": dict(k=6, mu=0.4791, se=0.1074, pct=61.5, ci=(22.5, 112.8), p=0.0066, tau2=0.0495, I2=75.4, Q=20.296),
    "SI": dict(k=6, mu=0.4207, se=0.0426, pct=52.3, ci=(36.5, 69.9), p=0.0002, tau2=0.0, I2=0.0, Q=4.823),
    "SEC": dict(k=5, mu=0.0362, se=0.1076, pct=3.7, ci=(-23.1, 39.8), p=0.7538, tau2=0.0538, I2=90.2, Q=40.625),
}

print("=" * 96)
print("第一步：复现 v4.0（必须与 meta_A.json 完全一致，否则停止）")
print("=" * 96)
ok_all = True
for q, rows in BASE.items():
    y = [r[1] for r in rows]; v = [r[2] ** 2 for r in rows]
    got = summarize(dl_kh(y, v), q)
    e = EXPECT[q]
    checks = [
        ("k", got["k"], e["k"]), ("mu", got["mu"], e["mu"]), ("se", got["se"], e["se"]),
        ("pct", got["pct"], e["pct"]), ("CI", got["ci"], e["ci"]), ("p", got["p"], e["p"]),
        ("tau2", got["tau2"], e["tau2"]), ("I2", got["I2"], e["I2"]), ("Q", got["Q"], e["Q"]),
    ]
    bad = []
    for name, g, x in checks:
        if isinstance(x, tuple):
            if abs(g[0] - x[0]) > 0.15 or abs(g[1] - x[1]) > 0.15: bad.append(f"{name}: {g} != {x}")
        elif isinstance(x, float):
            if abs(g - x) > max(0.0015, abs(x) * 0.006): bad.append(f"{name}: {g} != {x}")
        else:
            if g != x: bad.append(f"{name}: {g} != {x}")
    print(f"\n[{q}] k={got['k']} 合并={got['pct']}% CI{got['ci']} p={got['p']} tau2={got['tau2']} I2={got['I2']} Q={got['Q']}")
    if bad:
        ok_all = False
        for b in bad: print("   ✗ 不一致:", b)
    else:
        print("   ✓ 与 meta_A.json 完全一致")

print("\n" + "=" * 96)
print("复现结果：", "全部通过 → 可以加入新研究" if ok_all else "存在不一致 → 实现有误，禁止继续")
print("=" * 96)
