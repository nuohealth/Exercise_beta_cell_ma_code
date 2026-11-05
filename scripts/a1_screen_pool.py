import os, re, json, time

D = json.load(open("analysis/a1_pool_raw.json"))
meta, ab, counts = D["meta"], D["abs"], D["counts"]

ALREADY = {
 "17763840","19592624","22344613","22751691","25099941","37127822","30363954","31969136",
 "41103344","42050504","42698534",
}

EXW = r"exercis|training|physical activ|walking|aerobic|resist|hiit|interval|sedentar|sitting|lifestyle|qigong|yoga|tai chi|dance|step"
SEC = r"insulin secretion|c-?peptide|disposition index|first-phase|second-phase|first phase|second phase|airg|acute insulin response|beta-cell function|β-cell function|insulin secretory|insulin release|proinsulin|insulinogenic|homa2?-?\s?b|homa-?%?β|isr\b|gsis|secretion rate"
SEN = r"insulin sensitivity|\bsi\b|\bsi\(|isi|matsuda|m value|glucose infusion rate|\bgir\b|homa2?-?\s?ir|homa-?ir|clamp|1/fasting insulin|quicki"
DYN = r"clamp|ivgtt|intravenous glucose|deconvolut|minimal model|insulin secretion rate|\bisr\b|botnia|fsivgtt|ogtt.{0,40}deconvolut|oral minimal model|\bomm\b"
ANI = r"\b(rats?|mice|mouse|murine|rodent|porcine|canine|dogs?|zebrafish|in vitro|cell line|islets?)\b"
T1D = r"type 1 diabetes|type 1 diabetic|\bt1d\b|gestational diabetes|\bgdm\b"
DRUGT = r"metformin|sitagliptin|pioglitazone|rosiglitazone|empagliflozin|dapagliflozin|liraglutide|semaglutide|tirzepatide|exenatide|insulin glargine|sulfonylurea|gliclazide|glimepiride|alogliptin|linagliptin|saxagliptin|vildagliptin|canagliflozin|ertugliflozin|dulaglutide|orlistat|acarbose|repaglinide|nateglinide|miglitol|mitiglinide|curcumin|probiotic|synbiotic|vitamin d|cholecalciferol|fenugreek|berberine|resveratrol|carnosine|magnesium|omega-3|fish oil|creatine|yogurt|saffron|cinnamon|ginger"
REV = r"systematic review|meta-analysis|scoping review|narrative review|umbrella review|review of"

rows = []
for pid, m in meta.items():
    a = ab.get(pid, "")
    ti = m["title"] or ""
    low = (ti + " " + a).lower()
    tlow = ti.lower()
    pts = m["pt"] or []
    isrev = any(x in pts for x in ("Review", "Meta-Analysis", "Systematic Review")) or re.search(REV, tlow)
    rows.append(dict(
        pmid=pid, year=m["year"], journal=m["journal"], title=ti, pmc=m["pmc"],
        ab=a, groups=sorted(D["meta"].get(pid, {}).get("g", [])) if False else None,
        ex=bool(re.search(EXW, low)), sec=bool(re.search(SEC, low)), sen=bool(re.search(SEN, low)),
        dyn=bool(re.search(DYN, low)), ani=bool(re.search(ANI, tlow)), t1d=bool(re.search(T1D, low)),
        drug=bool(re.search(DRUGT, tlow)), rev=bool(isrev),
        trial=bool(re.search(r"randomi[sz]|controlled trial|\btrial\b|crossover|cross-over", low)),
        already=pid in ALREADY,
    ))

def dur_weeks(a):
    """从摘要里抓时长，返回周数或 None"""
    cands = []
    for mm in re.finditer(r"(\d+(?:\.\d+)?)\s*(week|weeks|wk|wks|month|months|mo|day|days|year|years)\b", a, re.I):
        v = float(mm.group(1)); u = mm.group(2).lower()
        if u.startswith("week") or u.startswith("wk"): cands.append(v)
        elif u.startswith("month") or u == "mo": cands.append(v * 4.345)
        elif u.startswith("day"): cands.append(v / 7)
        elif u.startswith("year"): cands.append(v * 52.14)
    return max(cands) if cands else None

for r in rows:
    r["weeks"] = dur_weeks(r["ab"])

# 初筛
keep = []
for r in rows:
    if r["already"] or r["rev"] or r["ani"] or r["t1d"] or r["drug"]: continue
    if not (r["ex"] and r["sec"]): continue
    if not r["trial"]: continue
    keep.append(r)

keep.sort(key=lambda r: (-(2 if r["dyn"] else 0) - (1 if r["sen"] else 0), r["year"]), reverse=False)
keep.sort(key=lambda r: (-((2 if r["dyn"] else 0) + (1 if r["sen"] else 0)), -int(r["year"] or 0)))

print("=" * 112)
print(f"池 {len(rows)} → 初筛通过 {len(keep)} 条（已排除：综述/动物/T1D-GDM/药名在题/已纳入）")
print("=" * 112)
for r in keep:
    tag = ("★DYN" if r["dyn"] else "  sur") + ("+SEN" if r["sen"] else "    ")
    w = f"{r['weeks']:.0f}w" if r["weeks"] is not None else " ? "
    print(f"{tag} {r['pmid']} {r['year']} {w:>4} {r['journal'][:32]:32} PMC={r['pmc'] or '-':12} {r['title'][:88]}")

json.dump(keep, open("analysis/a1_shortlist_auto.json", "w"), ensure_ascii=False, indent=1)
print(f"\n落盘 analysis/a1_shortlist_auto.json  ({len(keep)} 条)")
