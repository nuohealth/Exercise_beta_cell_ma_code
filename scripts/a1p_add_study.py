"""
A1' 第二步：加入 15677490（Diabetes 2005, IGT 日裔美国人，24 个月，饮食+耐力 vs 饮食+拉伸）
- 该研究是直接的随机两组对比（两组都吃饮食，差异=运动）→ 运动臂可分离，符合 §2.3
- 三个量都有 mean±SE：Si 4.7±0.5 vs 3.3±0.3；AIRg 294.9±50.0 vs 305.4±30.0；DI 作者只报 p=0.7（无数值）
- DI 用分量乘积导出：ln(DI)=ln(Si)+ln(AIRg)，SE 用独立量传播 —— 公式与出处必须写明
- 新层：IGT（糖耐量受损）→ 第三层，既非 T2D 也非正常糖耐量
"""
import json, math

def dl_kh(y, v):
    k = len(y); w = [1.0/x for x in v]; sw = sum(w); sw2 = sum(x*x for x in w)
    mu_fe = sum(wi*yi for wi, yi in zip(w, y))/sw
    Q = sum(wi*(yi-mu_fe)**2 for wi, yi in zip(w, y))
    C = sw - sw2/sw
    tau2 = max(0.0, (Q-(k-1))/C) if C > 0 else 0.0
    wr = [1.0/(vi+tau2) for vi in v]; swr = sum(wr)
    mu = sum(wi*yi for wi, yi in zip(wr, y))/swr
    Qre = sum(wi*(yi-mu)**2 for wi, yi in zip(wr, y))
    var_hk = Qre/((k-1)*swr) if k > 1 else 1.0/swr
    I2 = max(0.0, (Q-(k-1))/Q)*100 if Q > 0 else 0.0
    return dict(k=k, mu_ln=mu, se=math.sqrt(var_hk), tau2=tau2, I2=I2, Q=Q, df=k-1,
                pct=(math.exp(mu)-1)*100)

T = {1:12.706,2:4.303,3:3.182,4:2.776,5:2.571,6:2.447,7:2.365,8:2.306,9:2.262,
     10:2.228,11:2.201,12:2.179,13:2.160,14:2.145,15:2.131,16:2.120,17:2.110,
     18:2.101,19:2.093,20:2.086}
def tc(df): return T.get(df, 1.96)

def p_from_t(t, df):
    x = df/(df+t*t); a, b = df/2.0, 0.5
    if x <= 0: return 1.0
    lb = math.lgamma(a)+math.lgamma(b)-math.lgamma(a+b)
    def bcf(aa, bb, xx):
        MAXIT, EPS, FPMIN = 200, 3e-16, 1e-300
        qab, qap, qam = aa+bb, aa+1, aa-1
        c, d = 1.0, 1-qab*xx/qap
        if abs(d) < FPMIN: d = FPMIN
        d = 1/d; h = d
        for m in range(1, MAXIT+1):
            m2 = 2*m
            aa2 = m*(bb-m)*xx/((qam+m2)*(aa+m2))
            d = 1+aa2*d
            if abs(d) < FPMIN: d = FPMIN
            c = 1+aa2/c
            if abs(c) < FPMIN: c = FPMIN
            d = 1/d; h *= d*c
            aa2 = -(aa+m)*(qab+m)*xx/((aa+m2)*(qap+m2))
            d = 1+aa2*d
            if abs(d) < FPMIN: d = FPMIN
            c = 1+aa2/c
            if abs(c) < FPMIN: c = FPMIN
            d = 1/d; de = d*c; h *= de
            if abs(de-1) < EPS: break
        return h
    bt = math.exp(-lb+a*math.log(x)+b*math.log(1-x))
    if x < (a+1)/(a+b+2): return bt*bcf(a, b, x)/a
    return 1-bt*bcf(b, a, 1-x)/b

def rep(d):
    df = d["df"]; tcv = tc(df)
    lo = math.exp(d["mu_ln"]-tcv*d["se"]); hi = math.exp(d["mu_ln"]+tcv*d["se"])
    p = p_from_t(abs(d["mu_ln"]/d["se"]), df)
    return dict(k=d["k"], pct=round((math.exp(d["mu_ln"])-1)*100, 1),
                ci=[round((lo-1)*100, 1), round((hi-1)*100, 1)], p=round(p, 4),
                tau2=round(d["tau2"], 4), I2=round(d["I2"], 1), Q=round(d["Q"], 3),
                mu=round(d["mu_ln"], 4), se=round(d["se"], 4))

# ============ 15677490 效应量提取（透明记录） ============
m1, s1, n1 = 4.7, 0.5, 30      # Step2 + endurance
m2, s2, n2 = 3.3, 0.3, 32      # Step1 + stretching
ln_Si = math.log(m1/m2)
se_Si = math.sqrt((s1/m1)**2 + (s2/m2)**2)   # delta method，两组独立

a1m, a1s = 294.9, 50.0
a2m, a2s = 305.4, 30.0
ln_Ag = math.log(a1m/a2m)
se_Ag = math.sqrt((a1s/a1m)**2 + (a2s/a2m)**2)

ln_DI = ln_Si + ln_Ag
se_DI = math.sqrt(se_Si**2 + se_Ag**2)       # 和的两项独立（同臂但 Si/AIRg 由同一模型估计，保守视为独立）

print("="*100)
print("15677490 效应量（Diabetes 2005;54:340-6；IGT；24 个月；n=30 vs 32）")
print("="*100)
print(f"  Si   : {m1}±{s1} vs {m2}±{s2}  → ln={ln_Si:+.4f}  SE={se_Si:.4f}  = {(math.exp(ln_Si)-1)*100:+.1f}%")
print(f"  AIRg : {a1m}±{a1s} vs {a2m}±{a2s} → ln={ln_Ag:+.4f}  SE={se_Ag:.4f}  = {(math.exp(ln_Ag)-1)*100:+.1f}%")
print(f"  DI   : = Si × AIRg 导出       → ln={ln_DI:+.4f}  SE={se_DI:.4f}  = {(math.exp(ln_DI)-1)*100:+.1f}%")
print(f"         （作者本人只报 DI 组间 p=0.7，未给数值 → 这是本研究的关键矛盾）")
print(f"  DI 的 CI: {((math.exp(ln_DI-1.96*se_DI)-1)*100):+.1f}% ~ {((math.exp(ln_DI+1.96*se_DI)-1)*100):+.1f}%  ← 极宽")

# ============ 原始 8 arms ============
BASE = {
 "DI": [("DOSE-EX/HED",0.8628899551470398,0.16085335950863347),("DOSE-EX/MED",0.7178397931503168,0.16274808236086336),
        ("Karstoft/IWT",0.5080216964332565,0.13116726835138387),("Slentz/LOWMOD",0.5026162977469957,0.13552476903859137),
        ("Slentz/LOWVIG",0.19637502284164818,0.108532266286991),("Slentz/HIGHVIG",0.22098672646679193,0.0957828570939168),
        ("Lee2012/AE",None,None),("Lee2012/RE",None,None)],
 "SI": [("DOSE-EX/HED",0.6043159668533296,0.1561707556940531),("DOSE-EX/MED",0.4054651081081644,0.1575640236381383),
        ("Karstoft/IWT",0.40413088509502776,0.09746328437917222),("Slentz/LOWMOD",0.542574322080571,0.09889769019129266),
        ("Slentz/LOWVIG",0.3011050927839216,0.08820599358380829),("Slentz/HIGHVIG",0.40082472855166235,0.09353263171956532)],
 "SEC": [("DOSE-EX/HED",0.2468600779315258,0.06360982235087835),("DOSE-EX/MED",0.3220834991691132,0.06434423301016023),
         ("Slentz/LOWMOD",-0.06597680233078876,0.08555408044437353),("Slentz/LOWVIG",-0.14536883215782767,0.09234497349850307),
         ("Slentz/HIGHVIG",-0.20908844189875272,0.08720418418205479)],
}

# Lee 2012 各臂 SE：按 scripts/m2_stratified.py 的口径重算
#   ln = log((base+delta)/base)；se = se_delta / post
LEE = {
    "DI": {"base": {"AE": 538.9, "RE": 545.5},
           "delta": {"AE": (105.7, 78.0), "RE": (37.4, 77.1)}},
    "SI": {"base": {"AE": 2.2, "RE": 2.9},
           "delta": {"AE": (0.4, 0.2), "RE": (0.8, 0.2)}},
    "SEC": {"base": {"AE": 292.7, "RE": 213.9},
            "delta": {"AE": (3.9, 22.4), "RE": (-38.2, 21.5)}},
}

def lee_arm(q, arm):
    b = LEE[q]["base"][arm]; d, se_d = LEE[q]["delta"][arm]
    post = b + d
    return (math.log(post / b), se_d / post)

NEW_DI  = (ln_DI, se_DI)
NEW_SI  = (ln_Si, se_Si)
NEW_SEC = (ln_Ag, se_Ag)

def pool(axis, extra=True, lee_on=True):
    rows = []
    for name, y, se in BASE[axis]:
        if y is None: continue
        rows.append((name, y, se))
    if lee_on:
        for arm in ("AE", "RE"):
            y, se = lee_arm(axis, arm)
            rows.append(("Lee2012/" + arm, y, se))
    if extra:
        rows.append((("15677490" if axis!="SEC" else "15677490(AIRg)"),
                     NEW_DI[0] if axis=="DI" else (NEW_SI[0] if axis=="SI" else NEW_SEC[0]),
                     NEW_DI[1] if axis=="DI" else (NEW_SI[1] if axis=="SI" else NEW_SEC[1])))
    y = [r[1] for r in rows]; v = [r[2]**2 for r in rows]
    return rep(dl_kh(y, v)), [r[0] for r in rows]

print("\n"+"="*100)
print("重跑三轴（含/不含 15677490）")
print("="*100)
out = {}
for axis, lab in [("DI","Disposition index"),("SI","Insulin sensitivity"),("SEC","Insulin secretion")]:
    a, arms_a = pool(axis, extra=False)
    b, arms_b = pool(axis, extra=True)
    out[axis] = {"base": a, "new": b, "arms_new": arms_b}
    print(f"\n[{lab}]")
    print(f"  原 (k={a['k']}): {a['pct']:+.1f}% (95% CI {a['ci'][0]:+.1f} ~ {a['ci'][1]:+.1f}), p={a['p']}, I2={a['I2']}%, tau2={a['tau2']}")
    print(f"  新 (k={b['k']}): {b['pct']:+.1f}% (95% CI {b['ci'][0]:+.1f} ~ {b['ci'][1]:+.1f}), p={b['p']}, I2={b['I2']}%, tau2={b['tau2']}")
    print(f"  arms: {arms_b}")

json.dump(out, open("analysis/meta_A_a1p.json","w"), indent=1, ensure_ascii=False)
print("\n→ analysis/meta_A_a1p.json")
