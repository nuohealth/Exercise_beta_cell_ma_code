#!/usr/bin/env python3
"""
DOSE-EX (PMID 37127822, Nat Metab 2023) 对数可加分解预验证
原理: DI = Si x Secretion  =>  dlnDI = dlnSi + dlnSecretion  (恒等式, 非拟合)
残差 = dlnDI - (dlnSi + dlnSec)  -> 应接近0; 偏离=构造不匹配/清除率/四舍五入
数据源: PMC10229430 Table 2 (晚期DI=高糖钳夹; ISR; ISI)
"""
import math
rows = [
    # label, DI%, Si%, ISR%
    ("DCON vs CON (纯饮食减重)", 58, 18, 33),
    ("MED vs CON  (中量运动)",   105, 50, 38),
    ("HED vs CON  (高量运动)",   137, 83, 28),
    ("HED vs DCON (等减重加量)", 50, 55, -4),
]
print(f"{'对比':28s} {'dlnDI':>7s} {'dlnSi':>7s} {'dlnSec':>7s} {'残差':>7s} {'残差%':>6s} {'敏感性占比':>9s}")
for lab, di, si, sec in rows:
    a, b, c = math.log1p(di/100), math.log1p(si/100), math.log1p(sec/100)
    res = a - (b + c)
    print(f"{lab:28s} {a:7.4f} {b:7.4f} {c:7.4f} {res:+7.4f} {100*res/a:+6.1f}% {100*b/a:8.0f}%")
