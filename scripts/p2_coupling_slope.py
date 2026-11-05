#!/usr/bin/env python3
"""
P2/P3 预实验：运动与减重干预下"分泌-敏感性耦联斜率 b"的估计
================================================================
理论：DI = 分泌 x 敏感性  ->  dln(DI) = dln(Sec) + dln(Si)
      双曲线法则（纯代偿）预测：dln(Sec) = -1 * dln(Si)  =>  DI 不变
      故：   dln(Sec) = dln(DI) + b * dln(Si)
             b = -1  => 干预只沿曲线滑动（DI 不变，纯代偿）
             b =  0  => 敏感性升而分泌不动（DI 升高，曲线外移 = set point 上移）

数据来源（逐篇 PMC 全文表格核对，仅用原文报告的均值/几何均值比）：
  1) Slentz 2009 Diabetes Care  PMID 19592624  PMC2752909  Table 2（IVGTT 最小模型）
  2) DOSE-EX            Nat Metab 2023 PMID 37127822  PMC10229430 Table 2（高糖钳夹，GMR）
  3) Villareal 2008     Obesity      PMID 18388888  PMC2799929（C 肽最小模型）
  4) Pedersen 2014      Diabetologia PMID 25099941（三阶段高糖钳夹）
  5) Braun 2017         JAP          PMID 28473613（钳夹 + proinsulin）
  6) Zhang 2023         iScience     PMID 37485354  PMC10362261（10 周，体重未变）
  7) Lim 2011           Diabetologia PMID 21656330（8 周 600kcal，阳性对照）
"""
import math, statistics as st

def ln(x):  # 相对变化% -> 对数变化
    return math.log1p(x / 100.0)

# ---------- 1. Slentz 2009：四臂，相对基线的算术均值变化 ----------
# (% change)  DI, AIRg, Si
slentz = {
    "CON  (对照)":            (-65 / 1562 * 100, 22 / 490.6 * 100, -0.35 / 3.5 * 100),
    "LOMI (低量中强度)":      (742 / 1270 * 100, -10 / 457.3 * 100, 1.7 / 3.1 * 100),
    "LOVI (低量高强度)":      (255 / 1533 * 100, -43 / 445.5 * 100, 0.8 / 3.7 * 100),
    "HOVI (高量高强度)":      (255 / 1305 * 100, -77 / 505.6 * 100, 1.1 / 3.2 * 100),
}

def slope(points):
    """最小二乘 dlnSec ~ dlnSi"""
    xs = [p[1] for p in points]; ys = [p[0] for p in points]
    mx, my = st.mean(xs), st.mean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return sxy / sxx if sxx else float('nan')

print("=" * 78)
print("P2-A  Slentz 2009 (PMID 19592624)  四臂分解  [IVGTT 最小模型]")
print("=" * 78)
print(f"{'臂':22s} {'dlnDI':>8s} {'dlnSi':>8s} {'dlnAIRg':>9s} {'残差':>8s} {'b臂内':>7s}")
pts = []
for k, (di, air, si) in slentz.items():
    a, b_, c = ln(di), ln(si), ln(air)
    res = a - (b_ + c)
    pts.append((c, b_))
    print(f"{k:22s} {a:8.3f} {b_:8.3f} {c:9.3f} {res:+8.3f}")
print(f"\n  四臂回归斜率 b = {slope(pts):.3f}   (双曲线法则预测 = -1.000)")
print(f"  残差区间: {min(ln(d)- (ln(s)+ln(a)) for d,a,s in slentz.values()):+.3f} ~ "
      f"{max(ln(d)-(ln(s)+ln(a)) for d,a,s in slentz.values()):+.3f}  (系统性为正 -> 第三隔间/算术均值偏倚)")

# ---------- 2. DOSE-EX：三臂 vs CON，几何均值比 ----------
dosex = {
    "DCON (纯饮食减重)": (58, 18, 33),
    "MED  (中量运动)":   (105, 50, 38),
    "HED  (高量运动)":   (137, 83, 28),
}
print("\n" + "=" * 78)
print("P2-B  DOSE-EX (PMID 37127822)  三臂 vs 对照  [高糖钳夹, GMR]")
print("=" * 78)
print(f"{'臂':22s} {'dlnDI':>8s} {'dlnSi':>8s} {'dlnISR':>8s} {'残差':>8s} {'敏感性占比':>9s}")
pts2 = []
for k, (di, si, isr) in dosex.items():
    a, b_, c = ln(di), ln(si), ln(isr)
    pts2.append((c, b_))
    print(f"{k:22s} {a:8.3f} {b_:8.3f} {c:8.3f} {a-(b_+c):+8.3f} {100*b_/a:8.0f}%")
print(f"\n  三臂回归斜率 b = {slope(pts2):.3f}   (双曲线法则预测 = -1.000)")
print("  等减重那一格（HED vs DCON）: DeltaDI +50%, DeltaSi +55%, DeltaISR -4%  -> 敏感性占比 108%")

# ---------- 3. 单点研究（分泌 vs 敏感性变化 -> 隐含 b） ----------
print("\n" + "=" * 78)
print("P2-C  单点研究：敏感性升、分泌变  ->  隐含 b")
print("=" * 78)
single = {
    "Villareal 2008 (18388888)": (0.0, 116.0, 100.0, "C肽最小模型; ISR 不变, Si 翻倍"),
    "Pedersen 2014 (25099941)":  (0.0, 49.8, 66.2, "三阶段钳夹; 三种刺激下分泌均不变"),
    "Braun 2017 (28473613)":     (0.0, 25.0, None, "钳夹; proinsulin 不变(运动单用)"),
    "Zhang 2023 (37485354)":     (38.0, -6.3, 38.0, "体重未变, Matsuda 反降, 晚相DI +38%"),
}
print(f"{'研究':28s} {'dlnSec':>7s} {'dlnSi':>7s} {'dlnDI':>7s} {'隐含b':>7s}  说明")
for k, (sec, si, di, note) in single.items():
    a, b_ = ln(sec), ln(si)
    bb = a / b_ if abs(b_) > 0.02 else float('nan')
    print(f"{k:28s} {a:7.3f} {b_:7.3f} {(ln(di) if di else float('nan')):7.3f} {bb:7.3f}  {note}")

# ---------- 4. P3 阳性对照：纯热量限制 ----------
print("\n" + "=" * 78)
print("P3  阳性对照：纯热量限制/减重（不走运动）")
print("=" * 78)
print("  Lim 2011 (21656330): 600 kcal/d x 1 周")
print(f"    第一相胰岛素反应 0.19 -> 0.46 nmol/min/m2 = +142% (dlnSec = {ln(142.0):.3f})")
print("    肝胰岛素敏感性（HGO 抑制率）43% -> 74%")
print("  Steven 2016 (27002059): VLCD 8 周 -> 应答者的第一相胰岛素反应'回归'")
print("  DOSE-EX DCON 臂: ISR +33% 且三个干预臂 ISR 增幅相近（不减重剂量依赖）")
print("  => 方向确认：减重轴驱动分泌能力，运动轴驱动敏感性。H1 得到独立人群交叉支持。")

# ---------- 5. 汇总 ----------
print("\n" + "=" * 78)
print("汇总：运动 vs 双曲线法则")
print("=" * 78)
allb = [slope(pts), slope(pts2), 0.0, 0.0, 0.0]
print(f"  运动类 b 估计: Slentz {slope(pts):+.2f} | DOSE-EX {slope(pts2):+.2f} | "
      f"Villareal 0.00 | Pedersen 0.00 | Braun ~0.00")
print(f"  中位数 b ~ {st.median(allb):+.2f}   双曲线法则 b = -1.00")
print(f"  含义：敏感性每升高 10%，纯代偿时分泌应下降约 9.5%；实测运动试验中分泌几乎不动（~0-2%）。")
print(f"        即：运动不是在双曲线上'滑动'，而是把曲线本身外推（set point 上移）。")
