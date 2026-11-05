"""
A1' 第三步：分层 + 敏感性 + 交互检验（含 15677490）
关键判断：15677490 是 IGT（糖耐量受损）→ 既非 T2D 也非正常糖耐量，构成第三层
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
    return dict(k=k, mu_ln=mu, se=math.sqrt(var_hk), tau2=tau2, I2=I2, Q=Q, df=k-1)

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
            aa2 = m*(bb-m)*xx/((qam+m2)*(aa+m2)); d = 1+aa2*d
            if abs(d) < FPMIN: d = FPMIN
            c = 1+aa2/c
            if abs(c) < FPMIN: c = FPMIN
            d = 1/d; h *= d*c
            aa2 = -(aa+m)*(qab+m)*xx/((aa+m2)*(qap+m2)); d = 1+aa2*d
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
    p = p_from_t(abs(d["mu_ln"]/d["se"]), df) if d["se"] > 0 else float("nan")
    return dict(k=d["k"], pct=round((math.exp(d["mu_ln"])-1)*100, 1),
                ci=[round((lo-1)*100, 1), round((hi-1)*100, 1)], p=round(p, 4),
                tau2=round(d["tau2"], 3), I2=round(d["I2"], 1), Q=round(d["Q"], 2),
                mu=round(d["mu_ln"], 4), se=round(d["se"], 4), df=df)

# ---------- 臂表（study, arm, stratum, ln, se） ----------
ARMS = [
 # T2D 层
 ("DOSE-EX 2023","HED","t2d","DI",0.8628899551470398,0.16085335950863347),
 ("DOSE-EX 2023","HED","t2d","SI",0.6043159668533296,0.1561707556940531),
 ("DOSE-EX 2023","HED","t2d","SEC",0.2468600779315258,0.06360982235087835),
 ("DOSE-EX 2023","MED","t2d","DI",0.7178397931503168,0.16274808236086336),
 ("DOSE-EX 2023","MED","t2d","SI",0.4054651081081644,0.1575640236381383),
 ("DOSE-EX 2023","MED","t2d","SEC",0.3220834991691132,0.06434423301016023),
 ("Karstoft 2014","IWT","t2d","DI",0.5080216964332565,0.13116726835138387),
 ("Karstoft 2014","IWT","t2d","SI",0.40413088509502776,0.09746328437917222),
 # 正常糖耐量肥胖层
 ("Slentz 2009","LOWMOD","ngt","DI",0.5026162977469957,0.13552476903859137),
 ("Slentz 2009","LOWMOD","ngt","SI",0.542574322080571,0.09889769019129266),
 ("Slentz 2009","LOWMOD","ngt","SEC",-0.06597680233078876,0.08555408044437353),
 ("Slentz 2009","LOWVIG","ngt","DI",0.19637502284164818,0.108532266286991),
 ("Slentz 2009","LOWVIG","ngt","SI",0.3011050927839216,0.08820599358380829),
 ("Slentz 2009","LOWVIG","ngt","SEC",-0.14536883215782767,0.09234497349850307),
 ("Slentz 2009","HIGHVIG","ngt","DI",0.22098672646679193,0.0957828570939168),
 ("Slentz 2009","HIGHVIG","ngt","SI",0.40082472855166235,0.09353263171956532),
 ("Slentz 2009","HIGHVIG","ngt","SEC",-0.20908844189875272,0.08720418418205479),
]

LEE = {
 "DI": {"base":{"AE":538.9,"RE":545.5},"delta":{"AE":(105.7,78.0),"RE":(37.4,77.1)}},
 "SI": {"base":{"AE":2.2,"RE":2.9},"delta":{"AE":(0.4,0.2),"RE":(0.8,0.2)}},
 "SEC":{"base":{"AE":292.7,"RE":213.9},"delta":{"AE":(3.9,22.4),"RE":(-38.2,21.5)}},
}
for q, spec in LEE.items():
    for arm in ("AE","RE"):
        b = spec["base"][arm]; d, sed = spec["delta"][arm]; post = b+d
        ARMS.append(("Lee 2012", arm, "ngt", q, math.log(post/b), sed/post))

# 15677490（IGT 层）
LN_SI  = math.log(4.7/3.3);                 SE_SI  = math.sqrt((0.5/4.7)**2+(0.3/3.3)**2)
LN_AG  = math.log(294.9/305.4);             SE_AG  = math.sqrt((50.0/294.9)**2+(30.0/305.4)**2)
LN_DI  = LN_SI+LN_AG;                       SE_DI  = math.sqrt(SE_SI**2+SE_AG**2)
ARMS += [
 ("Kumagai/Shad 2005" if False else "15677490","endo","igt","DI",LN_DI,SE_DI),
 ("15677490","endo","igt","SI",LN_SI,SE_SI),
 ("15677490","endo","igt","SEC",LN_AG,SE_AG),
]

def pool(rows):
    y=[r[0] for r in rows]; v=[r[1]**2 for r in rows]
    return rep(dl_kh(y,v))

def sel(strata=None, q=None, excl_study=None, one_per_study=False):
    out=[]
    for s,arm,st,qq,ln,se in ARMS:
        if q and qq!=q: continue
        if strata and st not in strata: continue
        if excl_study and s in excl_study: continue
        out.append((s,ln,se))
    if one_per_study:
        seen={}
        for s,ln,se in out:
            if s not in seen or se < seen[s][1]: seen[s]=(ln,se)
        out=[(s,ln,se) for s,(ln,se) in seen.items()]
    return out

print("="*104)
print("分层结果（三轴 × 三层）")
print("="*104)
res={"T2D":{},"IGT":{},"NGT-obesity":{},"ALL":{}}
for q,lab in [("DI","Disposition index"),("SI","Insulin sensitivity"),("SEC","Insulin secretion")]:
    print(f"\n### {lab}")
    for st,name in [("t2d","T2D"),("igt","IGT"),("ngt","NGT-obesity"),(None,"ALL")]:
        rows=sel(strata=[st] if st else None, q=q)
        n = len(rows); nstud = len({s for s,_,_ in rows})
        if n==0: continue
        r = pool([(ln,se) for _,ln,se in rows])
        res[name][q]=dict(r, k_studies=nstud, arms=[s for s,_,_ in rows])
        if n==1:
            t=r["mu"]/r["se"]; p=p_from_t(abs(t), 1)
            print(f"  {name:12s} k=1 研究={nstud}  {r['pct']:+.1f}% (单一研究，无法合并)")
        else:
            print(f"  {name:12s} k={n:2d} 研究={nstud}  {r['pct']:+7.1f}%  CI[{r['ci'][0]:+7.1f},{r['ci'][1]:+7.1f}]  p={r['p']:.4f}  I²={r['I2']:5.1f}%  τ²={r['tau2']:.3f}")

print("\n"+"="*104)
print("敏感性分析")
print("="*104)
sens={}
for q in ("DI","SI","SEC"):
    print(f"\n### {q}")
    # 1) 全臂（主分析）
    r=pool([(ln,se) for _,ln,se in sel(q=q)]); sens[f"main_{q}"]=r
    print(f"  主分析              k={r['k']:2d}  {r['pct']:+7.1f}% CI[{r['ci'][0]:+7.1f},{r['ci'][1]:+7.1f}]")
    # 2) 一研究一臂
    rows=sel(q=q, one_per_study=True); r=pool([(ln,se) for _,ln,se in rows]); sens[f"one_per_study_{q}"]=r
    print(f"  一研究一臂          k={r['k']:2d}  {r['pct']:+7.1f}% CI[{r['ci'][0]:+7.1f},{r['ci'][1]:+7.1f}]  p={r['p']:.4f}")
    # 3) 剔除 15677490
    rows=sel(q=q, excl_study={"15677490"}); r=pool([(ln,se) for _,ln,se in rows]); sens[f"excl_15677490_{q}"]=r
    print(f"  剔除 15677490        k={r['k']:2d}  {r['pct']:+7.1f}% CI[{r['ci'][0]:+7.1f},{r['ci'][1]:+7.1f}]  p={r['p']:.4f}")
    # 4) 剔除 DOSE-EX
    rows=sel(q=q, excl_study={"DOSE-EX 2023"}); r=pool([(ln,se) for _,ln,se in rows])
    print(f"  剔除 DOSE-EX         k={r['k']:2d}  {r['pct']:+7.1f}% CI[{r['ci'][0]:+7.1f},{r['ci'][1]:+7.1f}]  p={r['p']:.4f}")

# ---------- 交互检验：T2D vs NGT（分泌轴），原稿 z=6.77 ----------
print("\n"+"="*104)
print("交互检验（分泌轴 T2D vs NGT），原稿报差异 +47.7%, z=6.77, p=0.001")
print("="*104)
sec_t2d = res["T2D"]["SEC"]; sec_ngt = res["NGT-obesity"]["SEC"]
diff = sec_t2d["mu"] - sec_ngt["mu"]
se_d = math.sqrt(sec_t2d["se"]**2 + sec_ngt["se"]**2)
z = diff/se_d
from math import erfc
pz = erfc(abs(z)/math.sqrt(2))
print(f"  T2D: {sec_t2d['pct']:+.1f}% ; NGT: {sec_ngt['pct']:+.1f}%")
print(f"  差异 = {(math.exp(diff)-1)*100:+.1f}%  (ln diff={diff:+.4f}, SE={se_d:.4f})")
print(f"  z = {z:.2f}   p = {pz:.4g}")
print(f"  （原稿: 差异 +47.7%, z=6.77, p=0.001 —— 含 IGT 第三层后两主层数值不变，故交互检验应一致）")

json.dump({"strata":res,"sens":sens,"interaction_Z":z,"interaction_p":pz,
           "diff_pct":(math.exp(diff)-1)*100},
          open("analysis/a1p_strata_sens.json","w"), indent=1, ensure_ascii=False)
print("\n→ analysis/a1p_strata_sens.json")
