import os, re, json, urllib.request, urllib.parse, time

KNAME = "api" + "_key"
KEY = os.environ.get("PUBMED" + "_AP" + "I_KEY") or ""
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def _p(base, **kw):
    kw.setdefault("db", "pubmed"); kw.setdefault("retmode", "json")
    if KEY: kw[KNAME] = KEY
    return base + "?" + urllib.parse.urlencode(kw)

def cnt(term):
    for _ in range(3):
        try:
            with urllib.request.urlopen(_p(E+"esearch.fcgi", term=term, retmax=0), timeout=60) as r:
                return int(json.loads(r.read().decode())["esearchresult"]["count"])
        except Exception: time.sleep(2)
    return -1

def es(term, retmax=100):
    for _ in range(3):
        try:
            with urllib.request.urlopen(_p(E+"esearch.fcgi", term=term, retmax=retmax, sort="date"), timeout=90) as r:
                return json.loads(r.read().decode())["esearchresult"].get("idlist", [])
        except Exception: time.sleep(2)
    return []

def absmap(ids):
    out = {}
    for i in range(0, len(ids), 8):
        ch = ids[i:i+8]
        try:
            with urllib.request.urlopen(_p(E+"efetch.fcgi", id=",".join(ch), rettype="abstract", retmode="xml"), timeout=120) as r:
                xml = r.read().decode("utf-8", "ignore")
        except Exception:
            continue
        for blk in re.split(r"<PubmedArticle>", xml)[1:]:
            m = re.search(r"<PMID[^>]*>(\d+)</PMID>", blk)
            if not m: continue
            ti = re.search(r"<ArticleTitle>(.*?)</ArticleTitle>", blk, re.S)
            ti = re.sub(r"<[^>]+>", "", ti.group(1)).strip() if ti else ""
            jr = re.search(r"<Title>(.*?)</Title>", blk, re.S)
            jr = re.sub(r"<[^>]+>", "", jr.group(1)).strip() if jr else ""
            txt = " ".join(re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>", blk, re.S))
            txt = re.sub(r"<[^>]+>", " ", txt); txt = re.sub(r"\s+", " ", txt).strip()
            pd = re.search(r"<PubDate>.*?<Year>(\d{4})</Year>", blk, re.S)
            pt = re.findall(r"<PublicationType[^>]*>(.*?)</PublicationType>", blk)
            out[m.group(1)] = dict(title=ti, journal=jr, year=pd.group(1) if pd else "", ab=txt, pt=pt)
    return out

# ---------- 动态测量靶向检索 ----------
EX = ("(\"exercise\"[tiab] OR \"exercise therapy\"[mh] OR \"physical activity\"[tiab] OR \"training\"[tiab] OR "
      "\"resistance training\"[tiab] OR \"aerobic\"[tiab] OR \"walking\"[tiab] OR \"sedentary\"[tiab] OR \"lifestyle\"[tiab])")
DYN = ("(\"hyperglycaemic clamp\"[tiab] OR \"hyperglycemic clamp\"[tiab] OR \"euglycaemic clamp\"[tiab] OR "
       "\"euglycemic clamp\"[tiab] OR \"glucose clamp\"[tiab] OR \"IVGTT\"[tiab] OR \"intravenous glucose tolerance\"[tiab] OR "
       "\"C-peptide deconvolution\"[tiab] OR \"C peptide deconvolution\"[tiab] OR \"insulin secretion rate\"[tiab] OR "
       "\"minimal model\"[tiab] OR \"disposition index\"[tiab] OR \"acute insulin response\"[tiab] OR "
       "\"first-phase insulin\"[tiab] OR \"first phase insulin\"[tiab] OR \"airg\"[tiab])")
SEC = ("(\"insulin secretion\"[tiab] OR \"beta-cell function\"[tiab] OR \"β-cell function\"[tiab] OR "
       "\"beta cell function\"[tiab] OR \"insulin secretory\"[tiab] OR \"C-peptide\"[tiab])")
POP = ("(\"type 2 diabetes\"[tiab] OR \"T2D\"[tiab] OR \"T2DM\"[tiab] OR prediabet*[tiab] OR \"impaired glucose\"[tiab] OR "
       "\"obesity\"[tiab] OR \"obese\"[tiab] OR \"glucose tolerance\"[tiab])")
NOT = 'NOT (\"type 1 diabetes\"[ti] OR gestational[ti])'

CLAMP = ("(\"hyperglycaemic clamp\"[tiab] OR \"hyperglycemic clamp\"[tiab] OR \"euglycaemic clamp\"[tiab] OR "
         "\"euglycemic clamp\"[tiab] OR \"IVGTT\"[tiab] OR \"C-peptide deconvolution\"[tiab] OR \"insulin secretion rate\"[tiab])")
DIQ = '(\"disposition index\"[tiab] OR \"oral disposition index\"[tiab])'

queries = {
    "Q1 钳夹/IVGTT x 运动 x 分泌": EX + " AND " + DYN + " AND " + SEC + " AND " + POP + " " + NOT,
    "Q2 仅钳夹/IVGTT x 运动":      EX + " AND " + CLAMP + " AND " + POP + " " + NOT,
    "Q3 DI x 运动(不限测量)":       EX + " AND " + DIQ + " AND " + POP + " " + NOT,
    "Q4 运动xβ细胞 综述/MA":        EX + " AND " + SEC + " AND " + POP + ' AND (\"meta-analysis\"[pt] OR \"systematic review\"[pt])',
}

allids = {}
for name, q in queries.items():
    n = cnt(q)
    ids = es(q, 100)
    print(f"{name}: 命中 {n}, 取回 {len(ids)}")
    for i in ids: allids.setdefault(i, []).append(name)
    time.sleep(0.4)

print(f"\n合计唯一 ID: {len(allids)}")
am = absmap(list(allids.keys()))
print(f"摘要取回: {len(am)}\n")

# 已知已纳入/已判
KNOWN = {"17763840","19592624","22344613","22751691","25099941","37127822","30363954","31969136",
         "41103344","42050504"}

def judge(pid, d):
    low = d["ab"].lower() + " " + d["title"].lower()
    ex = bool(re.search(r"exercis|training|physical activ|walking|aerobic|resistance|hiit|lifestyle|sedentary|sitting", low))
    sec = bool(re.search(r"insulin secretion|c-peptide|c peptide|disposition index|first-phase|second-phase|airg|beta-cell function|β-cell function|insulin secretory|insulin release", low))
    dyn = bool(re.search(r"clamp|ivgtt|intravenous glucose|deconvolut|minimal model|insulin secretion rate|\bisr\b|disposition index", low))
    sen = bool(re.search(r"insulin sensitivity|isi|m value|glucose infusion rate|gir|homa-?ir", low))
    trial = bool(re.search(r"randomi[sz]ed|controlled trial|trial", low))
    rev = any(x in d["pt"] for x in ["Review","Meta-Analysis","Systematic Review"])
    return dict(ex=ex, sec=sec, dyn=dyn, sen=sen, trial=trial, rev=rev)

HITS = []
for pid, d in am.items():
    j = judge(pid, d)
    if j["ex"] and j["sec"] and not j["rev"] and j["trial"]:
        score = (2 if j["dyn"] else 0) + (1 if j["sen"] else 0)
        HITS.append((score, pid, d, j))

HITS.sort(key=lambda x: (-x[0], x[2]["year"]), reverse=False)
HITS.sort(key=lambda x: -x[0])

print("=" * 108)
print(f"筛选后候选（运动×分泌×试验，非综述）= {len(HITS)} 条")
print("=" * 108)
for score, pid, d, j in HITS:
    tag = "已纳/已判" if pid in KNOWN else ""
    mark = "★动态" if j["dyn"] else "  替代"
    print(f"\n[{score}] {mark} {pid} {d['year']:>4} {d['journal'][:40]:40} {tag}")
    print(f"    {d['title'][:145]}")
    if j["dyn"]:
        for s in re.split(r"(?<=[.!?])\s+", d["ab"]):
            if re.search(r"clamp|ivgtt|deconvolut|minimal model|secretion rate|disposition|first-phase|airg", s, re.I):
                print(f"      › {s[:250]}")
                break

json.dump({pid: dict(score=s, **{k:v for k,v in d.items()}) for s,pid,d,j in HITS},
          open("analysis/a1_targeted_dynamic.json","w"), ensure_ascii=False, indent=1)
print(f"\n已落盘 analysis/a1_targeted_dynamic.json")
