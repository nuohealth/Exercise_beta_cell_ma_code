#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate supplementary materials S2-S5 and S8 from the frozen screening/extraction data.

Inputs (single source of truth):
  screening/w1_screen_scored.csv       868 records, W1 title/abstract verdicts
  screening/w2_L2_core.csv             113 records, secretion-outcome core set
  screening/w2_L2_fulltext.csv         266 records, L2 rule screen
  screening/w2_manual_verdicts.csv      42 records, manual full-text verdicts
  screening/w2_FINAL_included.csv        8 records, final include list
  screening/w2_FINAL_PDF_details.json    2 records recovered via PDF
  osf_upload/data/                     extraction table, decomposition, meta outputs
  figures/                             Fig1_PRISMA.png (300 dpi)

Outputs: supplement/
  S2_search_strings.md
  S3_screening_decisions.md
  S4_fulltext_verdicts.md
  S5_no_fulltext.md
  S8_data_and_code.md

Run: python3 scripts/build_supplement.py
"""
import csv
import glob
import json
import os
import shutil
import subprocess

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(PROJ, "supplement")
os.makedirs(OUT, exist_ok=True)


def rd(name):
    with open(os.path.join(PROJ, "screening", name), encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def esc(s):
    return (s or "").replace("|", "\\|").replace("\n", " ").strip()


# ---------------------------------------------------------------- S2
def build_s2():
    src = os.path.join(PROJ, "osf_upload", "Supplementary_Methods_S1_search_strings.md")
    txt = open(src, encoding="utf-8").read()
    txt = txt.replace(
        "# Supplementary Methods S1 — Search strings (as executed)",
        "# Supplementary Methods S2 — Search strings (as executed)",
    )
    txt = txt.replace("Recorded in Supplementary Table S6 (§ deviations D1).",
                      "Recorded in Supplementary Table S6 (item D2).")
    txt = txt.replace("(§ deviations D1)", "(item D2)")
    out = os.path.join(OUT, "S2_search_strings.md")
    open(out, "w", encoding="utf-8").write(txt)
    return out


# ---------------------------------------------------------------- S3
def build_s3():
    rows = rd("w1_screen_scored.csv")
    from collections import Counter
    vc = Counter(r["verdict"] for r in rows)
    ec = Counter(r["excl"] for r in rows if r["excl"])

    L = []
    L.append("# Supplementary Table S3 — Screening decisions for all records")
    L.append("")
    L.append("**Manuscript.** Exercise and pancreatic β-cell function in type 2 diabetes "
             "(decomposition meta-analysis)")
    L.append("")
    L.append("All unique records retrieved by the search, with the first-stage (title/abstract) "
             "verdict and, where applicable, the exclusion reason. Records marked "
             "`INCLUDE_SCREEN` proceeded to rule-based screening (Supplementary Table S4 chain).")
    L.append("")
    L.append(f"**Total records: {len(rows)}**")
    L.append("")
    L.append("## Verdict summary")
    L.append("")
    L.append("| Verdict | n | Meaning |")
    L.append("|---|---|---|")
    meanings = {
        "INCLUDE_SCREEN": "Passed title/abstract; proceeded to rule-based screening",
        "EXCLUDE_TOPIC": "Off topic at title/abstract",
        "EXCLUDE": "Excluded at title/abstract for the reason in the `excl` column",
        "NO_ABSTRACT": "No abstract available; could not be screened on abstract",
        "MAYBE_EX_PARTIAL": "Ambiguous or partial outcome; flagged for triage",
    }
    for k, v in vc.most_common():
        L.append(f"| `{k}` | {v} | {meanings.get(k, '')} |")
    L.append("")
    if ec:
        L.append("## Exclusion reasons (where `excl` is set)")
        L.append("")
        L.append("| Reason | n |")
        L.append("|---|---|")
        for k, v in ec.most_common():
            L.append(f"| {k} | {v} |")
        L.append("")

    L.append("## Record-level table")
    L.append("")
    L.append("| # | PMID | Year | Journal | Measure | Design | Verdict | Excl. reason | Title |")
    L.append("|---:|---|---:|---|---|---|---|---|---|")
    for i, r in enumerate(rows, 1):
        L.append(
            f"| {i} | {esc(r['pmid'])} | {esc(r['year'])} | {esc(r['journal'])} | "
            f"{esc(r['measure'])} | {esc(r['design'])} | `{esc(r['verdict'])}` | "
            f"{esc(r['excl']) or '—'} | {esc(r['title'])} |"
        )
    out = os.path.join(OUT, "S3_screening_decisions.md")
    open(out, "w", encoding="utf-8").write("\n".join(L) + "\n")
    return out, len(rows)


# ---------------------------------------------------------------- S4
def build_s4():
    core = rd("w2_L2_core.csv")
    verdicts = {r["pmid"]: r for r in rd("w2_manual_verdicts.csv")}
    included = {r["pmid"]: r for r in rd("w2_FINAL_included.csv")}
    xml = {os.path.basename(f).replace(".xml", "") for f in
           glob.glob(os.path.join(PROJ, "extraction", "fulltext", "*.xml"))}
    pdf2 = set(json.load(open(os.path.join(PROJ, "screening",
                                           "w2_FINAL_PDF_details.json"))).keys())

    read = [r for r in core if r["pmcid"] in xml or r["pmid"] in pdf2]
    noread = [r for r in core if r not in read]

    from collections import Counter
    vc = Counter(v["decision"] for v in verdicts.values())

    L = []
    L.append("# Supplementary Table S4 — Full-text verdicts")
    L.append("")
    L.append("**Manuscript.** Exercise and pancreatic β-cell function in type 2 diabetes "
             "(decomposition meta-analysis)")
    L.append("")
    L.append("The 113 records that reported any insulin-secretion outcome (the core set). "
             "Full text was retrievable for 44; each of those received an explicit verdict. "
             "The remaining 69 are listed in Supplementary Table S5.")
    L.append("")
    L.append(f"**Core set: {len(core)} · full text read: {len(read)} · not retrievable: {len(noread)}**")
    L.append("")
    L.append("## Verdict summary (42 manual verdicts + 2 PDF-recovered records whose "
             "verdicts are carried in the include list)")
    L.append("")
    L.append("| Verdict | n | Meaning |")
    L.append("|---|---|---|")
    vm = {
        "INCLUDE": "Eligible; contributes to the synthesis",
        "EXCLUDE": "Excluded after full-text reading",
        "MAYBE_LIFESTYLE": "Lifestyle/multi-component intervention; no separable exercise arm",
        "CONSTRUCT_EXAMPLE": "Meets inclusion criteria but its secretion outcomes are "
                             "algebraic functions of fasting covariates; retained outside "
                             "pooling as the H2 exemplar (§3.6)",
    }
    for k, v in vc.most_common():
        L.append(f"| `{k}` | {v} | {vm.get(k, '')} |")
    L.append("")

    L.append("## A. Studies included in the synthesis")
    L.append("")
    L.append("| # | PMID | Year | Journal | Weeks | Secretion construct | Clearance reported | Title |")
    L.append("|---:|---|---:|---|---:|---|---|---|")
    for i, r in enumerate(sorted(included.values(), key=lambda x: x["year"]), 1):
        L.append(f"| {i} | {esc(r['pmid'])} | {esc(r['year'])} | {esc(r['journal'])} | "
                 f"{esc(r['weeks'])} | {esc(r['tier'])} | {esc(r['clearance'])} | "
                 f"{esc(r['title'])} |")
    L.append("")

    L.append("## B. Full-text verdicts, all 42 manually adjudicated records")
    L.append("")
    L.append("| # | Decision | PMID | Year | Journal | Weeks | Reason | Title |")
    L.append("|---:|---|---|---:|---|---:|---|---|")
    for i, r in enumerate(sorted(verdicts.values(),
                                 key=lambda x: (x["decision"], x["year"])), 1):
        L.append(f"| {i} | `{esc(r['decision'])}` | {esc(r['pmid'])} | {esc(r['year'])} | "
                 f"{esc(r['journal'])} | {esc(r['weeks'])} | {esc(r['reason'])} | "
                 f"{esc(r['title'])} |")
    L.append("")

    L.append("## C. Core-set records with no retrievable full text")
    L.append("")
    L.append("Listed in full in Supplementary Table S5; summarised here as a count because "
             "no verdict could be assigned.")
    L.append("")
    L.append(f"**{len(noread)} records.** Their abstract-level assessment is in S5.")

    out = os.path.join(OUT, "S4_fulltext_verdicts.md")
    open(out, "w", encoding="utf-8").write("\n".join(L) + "\n")
    return out, len(core), len(read), len(noread)


# ---------------------------------------------------------------- S5
def build_s5():
    core = rd("w2_L2_core.csv")
    xml = {os.path.basename(f).replace(".xml", "") for f in
           glob.glob(os.path.join(PROJ, "extraction", "fulltext", "*.xml"))}
    pdf2 = set(json.load(open(os.path.join(PROJ, "screening",
                                           "w2_FINAL_PDF_details.json"))).keys())
    noread = [r for r in core if not (r["pmcid"] in xml or r["pmid"] in pdf2)]

    scored = {r["pmid"]: r for r in rd("w1_screen_scored.csv")}

    L = []
    L.append("# Supplementary Table S5 — Records with no retrievable full text")
    L.append("")
    L.append("**Manuscript.** Exercise and pancreatic β-cell function in type 2 diabetes "
             "(decomposition meta-analysis)")
    L.append("")
    L.append(f"{len(noread)} of the 113 core-set records had no open-access full text in "
             "PubMed Central or Europe PMC and no recoverable PDF. No full-text verdict "
             "could be assigned. Their abstract-level information is reproduced below so "
             "that a reader can judge what was missed.")
    L.append("")
    L.append("**These records are a stated limitation.** Abstracts suggest most would have "
             "been excluded on design grounds, but we cannot exclude the possibility that "
             "some report relevant data. Their construct provenance is entirely unknown — "
             "the proportion of formula-derived secretion indices in this unread literature "
             "cannot be estimated.")
    L.append("")
    L.append("| # | PMID | Year | Journal | Secretion construct (abstract-level) | Design | Title |")
    L.append("|---:|---|---:|---|---|---|---|")
    for i, r in enumerate(sorted(noread, key=lambda x: x["year"] or "0"), 1):
        s = scored.get(r["pmid"], {})
        L.append(f"| {i} | {esc(r['pmid'])} | {esc(r['year'])} | {esc(r['journal'])} | "
                 f"{esc(r['construct'])} | {esc(s.get('design'))} | {esc(r['title'])} |")
    out = os.path.join(OUT, "S5_no_fulltext.md")
    open(out, "w", encoding="utf-8").write("\n".join(L) + "\n")
    return out, len(noread)


# ---------------------------------------------------------------- S8
def build_s8():
    data = os.path.join(PROJ, "osf_upload", "data")
    files = sorted(os.listdir(data))
    sizes = {f: os.path.getsize(os.path.join(data, f)) for f in files}

    L = []
    L.append("# Supplementary S8 — Data, extraction table, analysis outputs and code")
    L.append("")
    L.append("**Manuscript.** Exercise and pancreatic β-cell function in type 2 diabetes "
             "(decomposition meta-analysis)")
    L.append("")
    L.append("**Registration.** OSF Registries, DOI 10.17605/OSF.IO/UXP7Q "
             "(registered 22 October 2026; retrospective with respect to data extraction).")
    L.append("")
    L.append("Everything below is publicly available at <https://osf.io/uxp7q/>. "
             "The registration snapshot is frozen; the same files are mirrored here so that "
             "the submission package is self-contained.")
    L.append("")
    L.append("## File manifest")
    L.append("")
    L.append("| File | Bytes | Contents |")
    L.append("|---|---:|---|")
    desc = {
        "数据提取表v1_定型_20261021.md": "Extraction instrument (v1). Defines every field, "
                                          "including the construct-source classification that "
                                          "the audit rests on.",
        "analysis_scripts.zip": "All analysis code: screening, decomposition, meta-analysis, "
                                 "figure generation, provenance scan.",
        "decomposition.csv": "Arm-level decomposition of Δln(DI) into sensitivity, secretion "
                              "and residual components (Table 2).",
        "gap_summary.json": "Evidence-gap summary over the core set.",
        "meta_A.json": "Primary random-effects meta-analysis outputs (DerSimonian–Laird τ², "
                        "Knapp–Hartung intervals), by axis.",
        "meta_B.json": "Sensitivity meta-analysis outputs (leave-one-out; construct-stratified).",
        "provenance_scan_v2.csv": "Construct-provenance scan across readable full texts: "
                                   "verdict, coefficient pairs, covariate lists, assay mentions.",
        "w1_prisma_counts.json": "PRISMA counts at the title/abstract stage.",
        "w2_FINAL_included.csv": "Final included studies with verdict tier.",
        "w2_FINAL_prisma.json": "Full PRISMA cascade counts, identification to inclusion.",
    }
    for f in files:
        L.append(f"| `{f}` | {sizes[f]:,} | {desc.get(f, '')} |")
    L.append("")
    L.append("## Reproducibility notes")
    L.append("")
    L.append("- **Single data source.** All pooled estimates derive from `decomposition.csv` "
             "and the two meta output files; no figure or table value was typed by hand.")
    L.append("- **Scale.** Decomposition is performed on the natural-log scale: "
             "Δln(DI) = Δln(secretion) + Δln(sensitivity) + residual. Percent changes are "
             "exp(Δln) − 1.")
    L.append("- **Estimator.** DerSimonian–Laird τ² with Knapp–Hartung adjusted confidence "
             "intervals throughout. No meta-regression and no Egger's test were performed "
             "(no power at k ≤ 6); both are pre-specified as not performed and this is "
             "stated in the Methods rather than left implicit.")
    L.append("- **No imputation.** Three eligible studies reported the outcomes without an "
             "extractable effect size or confidence interval. They are reported "
             "qualitatively and were not imputed into any pooled estimate.")
    L.append("- **Provenance classification.** The construct-source classifier was calibrated "
             "against three studies with known ground truth before being applied "
             "(Supplementary Table S6, item D9).")
    L.append("- **Software.** Python 3.12; NumPy, pandas, SciPy, Matplotlib. Scripts and "
             "pinned versions are inside `analysis_scripts.zip`.")
    L.append("")
    L.append("## Protocol deviations")
    L.append("")
    L.append("Every deviation from the registered protocol and every post hoc addition is "
             "listed in Supplementary Table S6 (22 items). None altered the three "
             "pre-specified hypotheses, the eligibility criteria, the outcome definitions, "
             "or the analysis plan.")
    out = os.path.join(OUT, "S8_data_and_code.md")
    open(out, "w", encoding="utf-8").write("\n".join(L) + "\n")
    return out, len(files)


def main():
    print("S2:", build_s2())
    print("S3:", build_s3())
    print("S4:", build_s4())
    print("S5:", build_s5())
    print("S8:", build_s8())

    # copy the two already-written supplements in
    for src, dst in (
        ("osf_upload/Supplementary_Table_S6_protocol_deviations.md",
         "S6_protocol_deviations.md"),
        ("osf_upload/Supplementary_S7_PRISMA_2020_checklist.md",
         "S7_PRISMA_2020_checklist.md"),
        ("osf_upload/protocol_methods_v1.md", "S1_protocol.md"),
    ):
        s = os.path.join(PROJ, src)
        if os.path.exists(s):
            shutil.copy2(s, os.path.join(OUT, dst))
            print("copied:", dst)

    # PDF render for submission systems that want supplements as one document later
    print("\nfiles in supplement/:")
    for f in sorted(os.listdir(OUT)):
        print(f"  {f:45s} {os.path.getsize(os.path.join(OUT, f)):>9,} B")


if __name__ == "__main__":
    main()
