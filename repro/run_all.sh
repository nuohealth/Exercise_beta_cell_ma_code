#!/usr/bin/env bash
# Reproduce every numeric output reported in the manuscript, and verify each
# against the archived copy in analysis/.
#
# Usage:  bash repro/run_all.sh
# Exit:   0 = all checks passed, 1 = at least one mismatch
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

PY="${PYTHON:-python3}"
FAILED=0

say() { printf '\n\033[1m== %s ==\033[0m\n' "$1"; }

run() {
  local script="$1"
  printf '  %-34s' "$(basename "$script")"
  if $PY "$script" > "/tmp/repro_$(basename "$script" .py).log" 2>&1; then
    echo "ran"
  else
    echo "FAILED TO RUN"
    sed 's/^/      /' "/tmp/repro_$(basename "$script" .py).log" | tail -12
    FAILED=1
  fi
}

# ===========================================================================
# Stage 1 — arm-level decomposition table (the pre-D27 eight-arm table).
#
# m1_decomposition.py emits arm-level values on the published (from-baseline)
# scale; v37_apply_net.py then converts the Slentz arms to net-of-control
# (se_net = sqrt(se_arm^2 + se_con^2)). Both steps are needed: running m1 alone
# gives a table that differs from the archived file. The v36 backup in analysis/
# is the frozen input to v37 and doubles as its idempotency guard.
# ===========================================================================
say "Stage 1 — arm-level decomposition (v3.6/v3.7 lineage)"
run scripts/m1_decomposition.py
run scripts/v37_apply_net.py

# ===========================================================================
# Stage 2 — eight-arm meta-analysis and the interaction test.
#
# meta_A_v313.json is not the manuscript's primary pooled estimate (that is the
# nine-arm A1' chain at Stage 4), but it is the input to the stratum x secretion
# interaction test, whose Knapp-Hartung t value (df = 5) is what the manuscript
# reports as "interaction p = 0.001". The interaction involves only the two
# secretion strata; Carr 2005 is in the IGT stratum and enters neither, so the
# test is unchanged by the later addition of that trial.
# ===========================================================================
say "Stage 2 — eight-arm meta-analysis + interaction test"
run scripts/m2_v313.py
run scripts/m4_interaction_test.py

# ===========================================================================
# Stage 3 — A1' chain: the pooled estimates actually reported in the manuscript.
#
# Adds Carr 2005 (PMID 15677490) to give 9 DI arms / 9 SI arms / 8 Sec arms
# from five studies across three glucose strata.
# ===========================================================================
say "Stage 3 — A1' meta-analysis (nine arms, five studies)"
run scripts/a1p_add_study.py
run scripts/a1p_meta.py
run scripts/a1p_strata_sens.py

# ===========================================================================
# Stage 4 — frozen arm-level dataset, sensitivity analyses, printed table.
# ===========================================================================
say "Stage 4 — arm set, sensitivity, ESM Table 2"
run scripts/arms_v500.py
run scripts/s7_sens_v500.py
run scripts/m2_table2_decomposition.py

# ===========================================================================
# Verification. Numeric outputs are compared byte-for-byte against the committed
# copies: the pipeline is deterministic, so any drift means the code and the
# archived result have diverged and the archived number is no longer
# reproducible.
# ===========================================================================
say "Verification — regenerated outputs vs committed copies"

FILES=(
  analysis/decomposition.csv
  analysis/decomposition_v313.csv
  analysis/meta_A_v313.json
  analysis/interaction_test_v314.json
  analysis/meta_A_a1p.json
  analysis/a1p_strata_sens.json
  analysis/arms_v500.json
  analysis/s7_sensitivity_v500.json
  analysis/table2_decomposition.csv
)

for f in "${FILES[@]}"; do
  printf '  %-42s' "$f"
  if [ ! -f "$f" ]; then
    echo "MISSING"; FAILED=1; continue
  fi
  if git -C "$ROOT" rev-parse --git-dir >/dev/null 2>&1; then
    if git -C "$ROOT" diff --quiet -- "$f"; then echo "PASS"
    else echo "FAIL (differs from committed copy)"; FAILED=1; fi
  else
    want="$(grep -E "  \Q$f\E$" repro/expected_md5.txt | awk '{print $1}')"
    got="$(md5sum "$f" | awk '{print $1}')"
    if [ "$want" = "$got" ]; then echo "PASS"
    else echo "FAIL (md5 $got, expected $want)"; FAILED=1; fi
  fi
done

# ===========================================================================
# Spot-check the headline numbers. A file that still hashes to itself proves
# the code is unchanged, not that the numbers are the ones in the paper, so
# every value quoted in the text is asserted explicitly here.
# ===========================================================================
say "Spot-check — values quoted in the manuscript"

$PY - <<'PYEOF'
import json, csv, sys
bad = 0

def chk(label, got, want, tol=0.051):
    global bad
    ok = (abs(got - want) <= tol) if isinstance(want, float) else (got == want)
    print(f"  {label:52s} {got!s:>9s} (paper {want})  {'PASS' if ok else 'FAIL'}")
    if not ok:
        bad += 1

# --- Table 1: pooled effects, nine-arm A1' chain -------------------------
m = json.load(open("analysis/meta_A_a1p.json"))
exp = {"DI": (9, 46.6, 19.5, 79.9), "SI": (9, 42.3, 28.2, 57.9), "SEC": (8, 0.2, -15.5, 18.9)}
for q, (k, pct, lo, hi) in exp.items():
    n = m[q]["new"]
    chk(f"Table 1 {q} k", n["k"], k)
    chk(f"Table 1 {q} pooled %", n["pct"], pct)
    chk(f"Table 1 {q} CI lower", n["ci"][0], lo)
    chk(f"Table 1 {q} CI upper", n["ci"][1], hi)

# --- Strata (Table 1, by glucose status) --------------------------------
s = json.load(open("analysis/a1p_strata_sens.json"))["strata"]
chk("T2D secretion % (+32.8)", s["T2D"]["SEC"]["pct"], 32.8)
chk("T2D secretion p (0.084)", s["T2D"]["SEC"]["p"], 0.0838)
chk("NGT secretion % (-10.0)", s["NGT-obesity"]["SEC"]["pct"], -10.0)
chk("NGT secretion p (0.073)", s["NGT-obesity"]["SEC"]["p"], 0.0725)
chk("T2D DI % (+97.2)", s["T2D"]["DI"]["pct"], 97.2)
chk("NGT DI % (+25.5)", s["NGT-obesity"]["DI"]["pct"], 25.5)

# --- Interaction test: the Knapp-Hartung t value the paper reports ------
it = json.load(open("analysis/interaction_test_v314.json"))["M1_subgroup_difference"]
chk("Interaction difference % (+47.7)", it["ratio_pct"], 47.7)
chk("Interaction KH p (0.001)", round(it["p_kh_t"], 3), 0.001)
chk("Interaction df (5)", it["df_t"], 5)

# --- Interaction from the nine-arm chain (no IGT arm enters it) ---------
chk("Interaction Z, nine-arm chain (6.77)", round(json.load(open("analysis/a1p_strata_sens.json"))["interaction_Z"], 2), 6.77)

# --- ESM Table 2: ten printed rows -------------------------------------
rows = list(csv.DictReader(open("analysis/table2_decomposition.csv")))
chk("ESM Table 2 row count (10)", len(rows), 10)
t2 = {r["arm"]: r for r in rows}
for arm, di, si, sec in [
    ("DOSE-EX HED", 0.863, 0.604, 0.247),
    ("Slentz high/vig", 0.221, 0.401, -0.209),
    ("Lee aerobic", 0.179, 0.167, 0.013),
    ("Carr 2005", 0.319, 0.354, -0.035),
]:
    chk(f"ESM T2 {arm} dlnDI", float(t2[arm]["dln_DI"]), di)
    chk(f"ESM T2 {arm} dlnSi", float(t2[arm]["dln_Si"]), si)
    chk(f"ESM T2 {arm} dlnSec", float(t2[arm]["dln_Sec"]), sec)

sys.exit(1 if bad else 0)
PYEOF
[ $? -ne 0 ] && FAILED=1

# ===========================================================================
echo
if [ "$FAILED" -eq 0 ]; then
  echo "ALL CHECKS PASSED"
else
  echo "SOME CHECKS FAILED — see logs in /tmp/repro_*.log"
fi
exit "$FAILED"
