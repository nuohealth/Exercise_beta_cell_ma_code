import os, re, json, urllib.request, urllib.parse, time

KNAME = "api" + "_key"
KEY = os.environ.get("PUBMED" + "_AP" + "I_KEY") or ""
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def _p(base, **kw):
    kw.setdefault("db", "pubmed"); kw.setdefault("retmode", "json")
    if KEY: kw[KNAME] = KEY
    return base + "?" + urllib.parse.urlencode(kw)

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
            # PMC id
            pmc = ""
            for pm in re.findall(r"<ArticleId IdType=\"pmc\">(PMC\d+)</ArticleId>", blk):
                pmc = pm
            out[m.group(1)] = dict(title=ti, journal=jr, year=pd.group(1) if pd else "", ab=txt, pt=pt, pmc=pmc)
    return out

SHORT = ["28473613","32503004","25882384","38513367","26258597","30307821","29865281",
         "31033904","40527119","37306399","25781537","21293446","36760521","37362426",
         "33290248","23704674","28490514","22717286","27480182","30939905","32816096"]

am = absmap(SHORT)
for pid in SHORT:
    d = am.get(pid)
    if not d:
        print(f"\n[{pid}] ❌ 未取回"); continue
    ab = d["ab"]
    print("\n" + "=" * 108)
    print(f"[{pid}] {d['year']} | {d['journal'][:55]} | PMC={d['pmc'] or '-'}")
    print(f"  {d['title']}")
    print(f"  PT: {[x for x in d['pt'] if x not in ('Journal Article','Research Support, N.I.H., Extramural','Research Support, Non-U.S. Gov\'t','Research Support, U.S. Gov\'t, Non-P.H.S.','Research Support, N.I.H., Intramural')]}")
    print(f"  摘要全文：\n  {ab}")
