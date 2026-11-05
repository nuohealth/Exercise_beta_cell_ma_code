#!/usr/bin/env python3
"""m6_weight_shares.py — 分层在合并权重中的占比（v3.18 增补）

目的
----
手稿、S6 D24†/D26† 与投稿信都引用了"正常糖层占处置指数轴合并权重的 65.5%"。
此数字此前只存在于文中，无脚本可复现。本脚本从 analysis/meta_A_v313.json 的
per_arm_data 直接求和，产出可核查的权重表，写入 analysis/weight_shares_v318.json。

单一数据源
----------
analysis/meta_A_v313.json -> per_arm_data[]，每臂含 weight_pct（DerSimonian–Laird
随机效应权重，已归一化到该轴合计 100.0）。

口径说明
--------
- weight_pct 是随机效应模型下该臂在**该结局轴**上的权重，不是样本量占比，
  也不是固定效应权重。三轴各归一化一次，故每轴合计均为 100.0。
- 分层归并：stratum == "T2D" 归 2 型糖尿病层；其余归正常糖耐量层。
- 处置指数轴 8 臂 / 4 研究；分泌轴 7 臂 / 3 研究（Karstoft 无分泌估计值，
  S6/S7 已披露，其 0.000 仅作分解占位）。
"""

import json
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "analysis", "meta_A_v313.json")
OUT = os.path.join(ROOT, "analysis", "weight_shares_v318.json")

AXIS_NAME = {"DI": "disposition index", "SI": "insulin sensitivity", "SEC": "insulin secretion"}


def main() -> None:
    with open(SRC, encoding="utf-8") as fh:
        data = json.load(fh)

    arms = data["per_arm_data"]
    result = {"source": "analysis/meta_A_v313.json", "axes": {}}

    for q in ("DI", "SI", "SEC"):
        rows = [a for a in arms if a["q"] == q]
        if not rows:
            continue
        total = sum(a["weight_pct"] for a in rows)

        by_stratum = defaultdict(lambda: {"arms": 0, "weight": 0.0, "studies": set()})
        for a in rows:
            key = "T2D" if a["stratum"] == "T2D" else "normoglycaemic obesity"
            by_stratum[key]["arms"] += 1
            by_stratum[key]["weight"] += a["weight_pct"]
            by_stratum[key]["studies"].add(a["study"])

        strata = {}
        for key, v in by_stratum.items():
            strata[key] = {
                "arms": v["arms"],
                "studies": len(v["studies"]),
                "study_names": sorted(v["studies"]),
                "weight_pct": round(v["weight"], 2),
                "share_pct": round(v["weight"] / total * 100, 2),
            }

        result["axes"][q] = {
            "axis": AXIS_NAME[q],
            "arms_total": len(rows),
            "studies_total": len({a["study"] for a in rows}),
            "weight_total": round(total, 2),
            "strata": strata,
            "per_arm": [
                {
                    "study": a["study"],
                    "arm": a["arm"],
                    "stratum": a["stratum"],
                    "weight_pct": a["weight_pct"],
                }
                for a in sorted(rows, key=lambda x: -x["weight_pct"])
            ],
        }

    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(result, fh, ensure_ascii=False, indent=2)

    # 控制台报告 + 断言（口径写死，改动即报错）
    print(f"[源] {SRC}")
    for q in ("DI", "SI", "SEC"):
        if q not in result["axes"]:
            continue
        ax = result["axes"][q]
        print(f"\n=== {q} ({ax['axis']}) — {ax['arms_total']} 臂 / {ax['studies_total']} 研究 ===")
        for name, s in sorted(ax["strata"].items(), key=lambda kv: -kv[1]["share_pct"]):
            print(
                f"  {name:24s} {s['arms']} 臂 {s['studies']} 研究 "
                f"权重 {s['weight_pct']:6.2f} = {s['share_pct']:5.2f}%"
            )
        # 归一化自检
        tot_share = sum(s["share_pct"] for s in ax["strata"].values())
        assert abs(tot_share - 100.0) < 0.05, f"{q} 分层占比合计 {tot_share} != 100"

    # 手稿引用的三个数（v3.14 起，S6 D24† 与投稿信同源）
    expect = {"DI": 65.5, "SI": 73.0, "SEC": 69.5}
    for q, want in expect.items():
        got = result["axes"][q]["strata"]["normoglycaemic obesity"]["share_pct"]
        assert abs(got - want) < 0.15, f"{q} 正常糖层占比 {got} != 文中 {want}"
    print("\n[自检] 正常糖层占比 65.5 / 73.0 / 69.5 三数与手稿、S6、投稿信一致")
    print(f"[输出] {OUT}")


if __name__ == "__main__":
    main()
