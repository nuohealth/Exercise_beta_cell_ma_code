# Exercise and pancreatic β-cell function across the glycaemic spectrum

**Reproducible analysis code for a systematic review and meta-analysis of randomized
controlled trials.**

The disposition index (DI = insulin secretion × insulin sensitivity) rises after
exercise training, but the decomposition shows the rise is carried almost entirely
by insulin sensitivity, not by insulin secretion.

- Protocol / registration: OSF Registration DOI `10.17605/OSF.IO/UXP7Q`
  (<https://osf.io/uxp7q/>)
- Manuscript: *Diabetes, Obesity and Metabolism* (submitted)

---

## 1. What this repository contains

| Folder | Contents |
|---|---|
| `scripts/` | All analysis, figure and manuscript-build code (Python 3.11+) |
| `screening/` | PRISMA-stage screening records: titles/abstracts, verdicts, final include list, PRISMA counts |
| `extraction/` | Per-study data provenance audit (values traced to table/figure/page in the source article) |
| `analysis/` | Intermediate and final numeric outputs (`*.json`, `*.csv`) produced by the scripts |
| `figures/` | Final publication figures (PNG, 300 dpi) |
| `repro/` | One-command reproduction script and expected-value checks |

Nothing here redistributes publisher PDFs. Bibliographic records and the numeric
values we extracted are included; the source articles must be obtained from the
publisher or an institutional subscription.

## 2. Requirements

```
Python >= 3.11
numpy, scipy, matplotlib, python-docx
```

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

`python-docx` is only needed to rebuild the submission DOCX
(`scripts/build_submission_docx_v74.py`); every statistical result and figure
can be reproduced without it.

## 3. Reproduce everything

```bash
bash repro/run_all.sh
```

The script runs the pipeline in dependency order, then compares every regenerated
output against the archived copy in `analysis/` and prints `PASS`/`FAIL` per file.
A clean run ends with:

```
ALL CHECKS PASSED
```

### Pipeline order (why this order)

```
m1_decomposition.py     Δln(DI) = Δln(Sec) + Δln(Si) + residual, per study arm
      ↓
m2_meta.py / m2_v313.py  DerSimonian–Laird τ² with Knapp–Hartung CIs, 3 axes
      ↓
m2_stratified.py         same models, stratified by baseline glycaemic state
      ↓
m4_interaction_test.py   test that the two strata differ in secretion direction
      ↓
arms_v500.py             single frozen arm-level dataset for all figures
      ↓
m6_weight_shares.py      stratum share of pooled weight (traceability of 65.5%)
      ↓
m3_*.py                  figures
```

## 4. Methods, precisely

**Effect size.** All meta-analyses are performed on the natural-log scale.
For each arm, ln(ratio) and its standard error are derived from the reported
ratio and 95% CI: `SE = [ln(hi) − ln(lo)] / (2 × 1.959964)`. Outcomes are
expressed as percent change on back-transformation.

**Pooling.** Random-effects model, DerSimonian–Laird estimator of τ², with
Knapp–Hartung adjustment of the confidence interval. Knapp–Hartung is not
optional here: with k = 7–8 the unadjusted Wald interval is anticonservative.

**Decomposition.** Because DI = Si × Sec, the identity
`Δln(DI) = Δln(Sec) + Δln(Si)` holds exactly per arm. The residual is reported
and used as a construction check (deviation signals mismatched constructs,
differential clearance, or rounding).

**Two contrasts.**
- Contrast A (primary): exercise vs no-intervention control.
- Contrast B: exercise added to matched weight loss vs weight loss alone —
  asks whether the exercise effect is independent of weight loss.

**Why k differs across axes.** DI and Si are available for 8 arms; Sec is
available for 7. One trial (Karstoft 2014) reports DI and Si but no usable
secretion effect size, so it contributes to two axes only. This is a data
availability fact, not an analytic choice — `arms_v500.py` encodes it once and
every figure reads from that file so the counts cannot drift between text
and display.

## 5. Headline results (regenerate with `repro/run_all.sh`)

| Axis | k | Δ | 95% CI | p | I² |
|---|---|---|---|---|---|
| Disposition index | 8 | +47.5% | 16.9 to 86.1 | 0.0055 | 75.3 |
| Insulin sensitivity | 8 | +42.7% | 26.7 to 60.7 | 0.0002 | 56.7 |
| Insulin secretion | 7 | +0.5% | −17.6 to 22.6 | 0.954 | 87.1 |

The secretion interval straddles zero with a point estimate of essentially no
change, while the sensitivity interval excludes zero by a wide margin. That
asymmetry is the paper's core claim.

Stratified analysis (`analysis/stratified.json`) shows the two strata disagree in
secretion direction: T2D +32.8% versus normoglycaemic obesity −10.0%. The
interaction test (`analysis/interaction_test_v314.json`) quantifies how far that
opposition is from chance; see the manuscript for the pre-specified interpretation
limits.

## 6. Data sources

Studies were identified by systematic search of PubMed, Europe PMC and reference
lists; search strings are archived in the OSF registration and in
`screening/`. Every extracted value carries a provenance record naming the
table, figure or page it came from (`extraction/provenance_scan_v2.csv`).

The screening pool files (`screening/w1_*.csv/json`) contain bibliographic
metadata retrieved from public APIs. They are included so that the PRISMA
counts can be re-derived without repeating the retrieval.

## 7. Known limitations of the code

- Some scripts import sibling modules (`m2_meta`, `m3_figures_v2`) and rely on
  `sys.path` being extended with the `scripts/` directory, which they do
  themselves. Invoke them as `python3 scripts/xxx.py` from the repository root;
  `repro/run_all.sh` already does this.
- Figures are drawn with matplotlib defaults plus explicit rcParams; exact
  pixel-level identity across matplotlib major versions is not guaranteed.
  Numeric outputs are version-stable and are the primary check.
- `build_submission_docx_v74.py` points at a manuscript path in the working
  directory; it is build tooling for our submission, not part of the
  statistical reproduction.

## 8. License

Code: MIT (`LICENSE`). Derived numeric data and figures: CC BY 4.0.
Third-party bibliographic records remain subject to their sources' terms.

## 9. Contact

Zhongwei Sun — sunzhongwei@gzgs.edu.cn
School of Physical Education, Guangzhou College of Technology and Business
