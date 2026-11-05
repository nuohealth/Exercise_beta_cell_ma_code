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
Stage 1  m1_decomposition.py    Δln(DI) = Δln(Sec) + Δln(Si) + residual, per arm
         v37_apply_net.py       rescale the Slentz arms to net-of-control
         ↓
Stage 2  m2_v313.py             eight-arm models (v313 lineage)
         m4_interaction_test.py stratum × secretion interaction (KH t, df = 5)
         ↓
Stage 3  a1p_add_study.py       add Carr 2005 → nine arms, three strata
         a1p_meta.py            pooled estimates reported in the paper
         a1p_strata_sens.py     strata, sensitivity, interaction
         ↓
Stage 4  arms_v500.py           single frozen arm set behind every figure
         s7_sens_v500.py        leave-one-out, one-arm-per-study, weight shares
         m2_table2_decomposition.py  the ten printed rows of ESM Table 2
         ↓
         m3_*.py                figures
```

Two lineages are both present on purpose. The **A1′ chain (Stage 3)** produces the
pooled estimates reported in the manuscript. The **v313 chain (Stage 2)** is
retained because its eight-arm model is the input to the interaction test quoted
in the text; the test involves only the two secretion strata, and Carr 2005 sits
in the impaired-glucose-tolerance stratum, so adding that trial leaves the test
unchanged. Running only one lineage would leave part of the paper unreproducible.

## 4. Methods, precisely

**Effect size.** All meta-analyses are performed on the natural-log scale.
For each arm, ln(ratio) and its standard error are derived from the reported
ratio and 95% CI: `SE = [ln(hi) − ln(lo)] / (2 × 1.959964)`. Outcomes are
expressed as percent change on back-transformation.

**Pooling.** Random-effects model, DerSimonian–Laird estimator of τ², with
Knapp–Hartung adjustment of the confidence interval. Knapp–Hartung is not
optional here: with k = 8–9 the unadjusted Wald interval is anticonservative.

**Decomposition.** Because DI = Si × Sec, the identity
`Δln(DI) = Δln(Sec) + Δln(Si)` holds exactly per arm. The residual is reported
and used as a construction check (deviation signals mismatched constructs,
differential clearance, or rounding).

**Two contrasts.**
- Contrast A (primary): exercise vs no-intervention control.
- Contrast B: exercise added to matched weight loss vs weight loss alone —
  asks whether the exercise effect is independent of weight loss.

**Why k differs across axes.** DI and Si are available for 9 arms; Sec is
available for 8. One trial (Karstoft 2014) reports DI and Si but no usable
secretion effect size, so it contributes to two axes only. This is a data
availability fact, not an analytic choice — `arms_v500.py` encodes it once and
every figure reads from that file so the counts cannot drift between text
and display.

**Shared controls.** Four DOSE-EX rows, three Slentz rows and two Lee rows come
from trials whose arms share a single control. Each arm's variance is inflated by
the shared-control structure (Slentz: se × √2, carried in the standard errors),
and the diabetes stratum's variance is understated in the usual way. The
manuscript states this and reports the interaction as exploratory for that
reason.

## 5. Headline results (regenerate with `repro/run_all.sh`)

Nine arms from five studies; eight arms carry a secretion estimate, because
Karstoft 2014 reports no usable secretion effect size.

| Axis | k | Δ | 95% CI | p | I² |
|---|---|---|---|---|---|
| Disposition index | 9 | +46.6% | +19.5 to +79.9 | 0.0026 | 71.8 |
| Insulin sensitivity | 9 | +42.3% | +28.2 to +57.9 | 0.0001 | 50.6 |
| Insulin secretion | 8 | +0.2% | −15.5 to +18.9 | 0.978 | 85.1 |

The secretion interval straddles zero with a point estimate of essentially no
change, while the sensitivity interval excludes zero by a wide margin. That
asymmetry is the paper's core claim — and it is deliberately *not* stated as
evidence that secretion is unaffected, because a near-null pooled over strata of
opposite sign averages them away.

**The near-null is an average of two opposite stratum effects.** By baseline
glucose status, secretion rose **+32.8% (−17.6 to +114.2, p = 0.084)** in the two
arms from type 2 diabetes and fell **−10.0% (−20.3 to +1.5, p = 0.073)** in the
five arms from obesity with normal glucose tolerance. The difference is +47.7%
(exploratory interaction, Knapp–Hartung t, df = 5, **p = 0.001**). Neither
stratum is individually significant; what the data support is the opposition,
not either stratum effect taken alone. Arms sharing one control understate
variance, so the contrast is reported as exploratory throughout.

The impaired-glucose-tolerance stratum rests on a single trial (Carr 2005) and is
reported without pooling.

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
- **Build tooling is included for transparency, not for execution.**
  `build_submission_docx_v74.py`, `build_supplement.py`, `check_v74_dom.py` and
  `count_words.py` generate the submission files and the supplement from our
  internal working tree. `build_supplement.py` reads
  `osf_upload/Supplementary_Methods_S1_search_strings.md` and
  `extraction/fulltext/*.xml`, neither of which is redistributed here — the
  former is unpublished protocol material and the latter is publisher
  full text. Running it will raise `FileNotFoundError`; that is expected and
  does not affect any statistic. Every number in the manuscript is reproduced
  by `repro/run_all.sh` without these scripts.
- **Two weight-share files, and only one is current.** `analysis/weight_shares_v318.json`
  holds the *superseded* eight-arm shares (normal-glucose stratum 65.5%). The
  figure the manuscript reports (67.8%) comes from `analysis/s7_sensitivity_v500.json`,
  computed on the nine-arm set by `scripts/s7_sens_v500.py`. The eight-arm file is
  kept because the manuscript's change log refers to it explicitly ("not the
  65.5% stated in the previous version"). If you are checking a number in the
  submitted paper, use the `s7_` file.
- `build_submission_docx_v74.py` points at a manuscript path in the working
  directory; it is build tooling for our submission, not part of the
  statistical reproduction.

## 8. License

Code: MIT (`LICENSE`). Derived numeric data and figures: CC BY 4.0.
Third-party bibliographic records remain subject to their sources' terms.

## 9. Contact

Zhongwei Sun — sunzhongwei@gzgs.edu.cn
School of Physical Education, Guangzhou College of Technology and Business
