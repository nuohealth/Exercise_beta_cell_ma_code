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

def es_all(term, cap=400):
    """分页取全部 ID"""
    ids, retstart = [], 0
    while retstart < cap:
        try:
            with urllib.request.urlopen(_p(E+"esearch.fcgi", term=term, retmax=200,
                                            retstart=retstart, sort="date"), timeout=90) as r:
                d = json.loads(r.read().decode())["esearchresult"]
        except Exception:
            time.sleep(2); continue
        batch = d.get("idlist", [])
        ids += batch
        total = int(d.get("count", 0))
        retstart += 200
        if not batch or retstart >= total: break
        time.sleep(0.35)
    return ids, total

def summ(ids):
    """esummary → 权威 PMCID + 期刊 + 年 + 类型"""
    out = {}
    for i in range(0, len(ids), 50):
        ch = ids[i:i+50]
        for attempt in range(3):
            try:
                with urllib.request.urlopen(_p(E+"esummary.fcgi", id=",".join(ch)), timeout=90) as r:
                    d = json.loads(r.read().decode())["result"]
                break
            except Exception:
                time.sleep(2); d = None
        if not d: continue
        for pid in ch:
            u = d.get(pid)
            if not u or pid == "uids": continue
            pmc = ""
            aids = u.get("articleids", [])
            for a in aids:
                if a.get("idtype") == "pmc":
                    pmc = a.get("value", "")
            out[pid] = dict(title=u.get("title", ""), journal=u.get("source", ""),
                            year=(u.get("pubdate", "") or "")[:4], pmc=pmc,
                            pt=u.get("pubtype", []))
        time.sleep(0.35)
    return out

def absmap(ids):
    out = {}
    for i in range(0, len(ids), 8):
        ch = ids[i:i+8]
        for attempt in range(3):
            try:
                with urllib.request.urlopen(_p(E+"efetch.fcgi", id=",".join(ch),
                                                rettype="abstract", retmode="xml"), timeout=120) as r:
                    xml = r.read().decode("utf-8", "ignore")
                break
            except Exception:
                time.sleep(3); xml = ""
        if not xml: continue
        for blk in re.split(r"<PubmedArticle>", xml)[1:]:
            m = re.search(r"<PMID[^>]*>(\d+)</PMID>", blk)
            if not m: continue
            txt = " ".join(re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>", blk, re.S))
            txt = re.sub(r"<[^>]+>", " ", txt); txt = re.sub(r"\s+", " ", txt).strip()
            out[m.group(1)] = txt
        time.sleep(0.3)
    return out

EX = ('("exercise"[tiab] OR "exercise therapy"[mh] OR "physical activity"[tiab] OR "training"[tiab] OR '
      '"resistance training"[tiab] OR "aerobic"[tiab] OR "walking"[tiab] OR "sedentary"[tiab] OR '
      '"lifestyle intervention"[tiab] OR "lifestyle"[tiab] OR "exercise training"[tiab])')
CLAMP = ('("hyperglycaemic clamp"[tiab] OR "hyperglycemic clamp"[tiab] OR "euglycaemic clamp"[tiab] OR '
         '"euglycemic clamp"[tiab] OR "glucose clamp"[tiab] OR "IVGTT"[tiab] OR "intravenous glucose tolerance"[tiab] OR '
         '"C-peptide deconvolution"[tiab] OR "C peptide deconvolution"[tiab] OR "insulin secretion rate"[tiab] OR '
         '"minimal model"[tiab] OR "Botnia clamp"[tiab] OR "botnia"[tiab])')
DIQ = '("disposition index"[tiab] OR "oral disposition index"[tiab] OR "adaptation index"[tiab])'
SEC = ('("insulin secretion"[tiab] OR "beta-cell function"[tiab] OR "beta cell function"[tiab] OR '
       '"insulin secretory"[tiab] OR "C-peptide"[tiab] OR "proinsulin"[tiab])')
POP = ('("type 2 diabetes"[tiab] OR "T2D"[tiab] OR "T2DM"[tiab] OR prediabet*[tiab] OR "impaired glucose"[tiab] OR '
       '"impaired fasting"[tiab] OR "obesity"[tiab] OR "obese"[tiab] OR "glucose tolerance"[tiab] OR '
       '"metabolic syndrome"[tiab] OR "overweight"[tiab])')
NOTX = 'NOT ("type 1 diabetes"[ti] OR "type 1"[ti] OR gestational[ti] OR "animal"[ti] OR mice[ti] OR rat[ti])'

QS = {
    "A 钳夹/IVGTT×运动×人群": EX + " AND " + CLAMP + " AND " + POP + " " + NOTX,
    "B DI×运动×人群":          EX + " AND " + DIQ + " AND " + POP + " " + NOTX,
    "C 运动×分泌×人群":        EX + " AND " + SEC + " AND " + POP + " " + NOTX,
}

pool = {}
counts = {}
for name, q in QS.items():
    ids, total = es_all(q)
    counts[name] = (total, len(ids))
    print(f"{name}: 命中 {total}, 取回 {len(ids)}", flush=True)
    for i in ids:
        pool.setdefault(i, set()).add(name[0])

print(f"\n合并唯一: {len(pool)}", flush=True)
meta = summ(list(pool.keys()))
print(f"esummary 取回: {len(meta)}", flush=True)
ab = absmap(list(meta.keys()))
print(f"摘要取回: {len(ab)}\n", flush=True)

json.dump({"counts": counts, "meta": meta, "abs": ab},
          open("analysis/a1_pool_raw.json", "w"), ensure_ascii=False, indent=1)
print("落盘 analysis/a1_pool_raw.json")
