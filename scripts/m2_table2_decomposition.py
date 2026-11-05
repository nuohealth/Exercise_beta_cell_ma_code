#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regenerate ESM Table 2 (arm-level decomposition of the change in disposition index).

Single data source: analysis/arms_v500.json (nine DI arms, nine SI arms, eight Sec arms).
Output:             analysis/table2_decomposition.csv

Why this script exists
----------------------
The archived v4.0 file `analysis/decomposition.csv` holds the **pre-D27, eight-arm**
decomposition (four studies, no Lee 2012, no Carr 2005, includes Eriksen 2007 and the
Slentz reference arm). It is retained unchanged for provenance but it is NOT the table
printed as ESM Table 2 in the submitted manuscript: the submitted table has ten rows —
the nine pooled arms plus the DOSE-EX diet-only reference arm — and includes Lee 2012
and Carr 2005.

ESM Table 2 is therefore reproduced here from the same frozen arm set the pooled
estimates and figures use, so that the printed table and the pooled numbers cannot
drift apart.

Conventions carried from the manuscript
---------------------------------------
  Δln(DI) = Δln(Sec) + Δln(Si) + residual, per arm, on the natural-log scale.
  Shares are ratios to Δln DI and are undefined (left blank) where Δln DI is
  mechanical or near zero:
    - Carr 2005: DI is reconstructed as its own product, so residual is 0 by
      construction and the shares would imply an informative decomposition.
    - Karstoft 2014: no secretion estimate (NE), so Sec share is NE.
  The DOSE-EX diet-only arm (DCON) is a reference row, not an effect size.
"""
import csv
import json
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(BASE, "analysis", "arms_v500.json")
OUT = os.path.join(BASE, "analysis", "table2_decomposition.csv")

# Display order and labels exactly as printed in ESM Table 2.
ROWS = [
    ("DOSE-EX 2023", "DCON",    "DOSE-EX DCON",     "T2D"),
    ("DOSE-EX 2023", "MED",     "DOSE-EX MED",      "T2D"),
    ("DOSE-EX 2023", "HED",     "DOSE-EX HED",      "T2D"),
    ("Karstoft 2014", "IWT",    "Karstoft IWT",     "T2D"),
    ("Slentz 2009", "LOWMOD",   "Slentz low/mod",   "NGT"),
    ("Slentz 2009", "LOWVIG",   "Slentz low/vig",   "NGT"),
    ("Slentz 2009", "HIGHVIG",  "Slentz high/vig",  "NGT"),
    ("Lee 2012", "AE",          "Lee aerobic",      "NGT"),
    ("Lee 2012", "RE",          "Lee resistance",   "NGT"),
    ("Carr 2005", "",           "Carr 2005",        "IGT"),
]

# Rows whose decomposition is mechanical rather than informative.
MECHANICAL = {("Carr 2005", "")}

# The diet-only reference arm: present in the table, absent from every pooled arm set.
REFERENCE = {("DOSE-EX 2023", "DCON")}


def load():
    arms = json.load(open(SRC, encoding="utf-8"))
    idx = {}
    for q in ("DI", "SI", "SEC"):
        for r in arms[q]:
            idx[(r["study"], r["arm"], q)] = r["ln"]
    # The diet-only arm enters no pooled analysis, so it is absent from arms_v500.json.
    # Read it from the archived arm-level table rather than typing the values here.
    dcon = {}
    arch = os.path.join(BASE, "analysis", "decomposition.csv")
    for row in csv.DictReader(open(arch, encoding="utf-8")):
        if row["study"] == "DOSE-EX 2023" and row["arm"] == "DCON":
            dcon = {"DI": float(row["dln_DI"]),
                    "SI": float(row["dln_Si"]),
                    "SEC": float(row["dln_Sec"])}
            break
    if not dcon:
        raise SystemExit("DCON reference row not found in analysis/decomposition.csv")
    return idx, dcon


def main():
    idx, dcon = load()
    out_rows = []

    for study, arm, label, stratum in ROWS:
        if (study, arm) in REFERENCE:
            di, si, sec = dcon["DI"], dcon["SI"], dcon["SEC"]
        else:
            di = idx.get((study, arm, "DI"))
            si = idx.get((study, arm, "SI"))
            sec = idx.get((study, arm, "SEC"))

        mechanical = (study, arm) in MECHANICAL
        # Karstoft contributes no secretion estimate: residual is the gap to Si alone.
        sec_missing = sec is None
        resid = None
        if di is not None and si is not None:
            if sec_missing:
                resid = di - si
            else:
                resid = di - si - sec

        if mechanical or di is None or si is None or abs(di) < 1e-9:
            sh_si = sh_sec = None
        else:
            sh_si = si / di
            sh_sec = None if sec_missing else sec / di

        out_rows.append(dict(
            arm=label, stratum=stratum,
            dln_DI=di, dln_Si=si, dln_Sec=sec, residual=resid,
            share_Si=sh_si, share_Sec=sh_sec,
            reference="yes" if (study, arm) in REFERENCE else "",
            mechanical="yes" if mechanical else "",
        ))

    def f(x, nd=3):
        return "" if x is None else f"{x:+.{nd}f}"

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["arm", "stratum", "dln_DI", "dln_Si", "dln_Sec",
                    "residual", "share_Si", "share_Sec", "reference", "mechanical"])
        for r in out_rows:
            w.writerow([r["arm"], r["stratum"], f(r["dln_DI"]), f(r["dln_Si"]),
                        "" if r["dln_Sec"] is None else f(r["dln_Sec"], 3),
                        "" if r["residual"] is None else f(r["residual"], 3),
                        "" if r["share_Si"] is None else f"{r['share_Si']:.2f}",
                        "" if r["share_Sec"] is None else f"{r['share_Sec']:.2f}",
                        r["reference"], r["mechanical"]])

    print(f"{'arm':18s}{'stratum':9s}{'dDI':>9s}{'dSi':>9s}{'dSec':>9s}{'resid':>9s}  shares")
    for r in out_rows:
        print(f"{r['arm']:18s}{r['stratum']:9s}{f(r['dln_DI']):>9s}{f(r['dln_Si']):>9s}"
              f"{f(r['dln_Sec']):>9s}{f(r['residual']):>9s}  "
              f"{'' if r['share_Si'] is None else round(r['share_Si'],2)} / "
              f"{'' if r['share_Sec'] is None else round(r['share_Sec'],2)}")
    print(f"\n→ {os.path.relpath(OUT, BASE)}")


if __name__ == "__main__":
    main()
