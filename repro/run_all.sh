#!/usr/bin/env bash
# Reproduce every numeric output in this repository and verify it against the
# archived copies in analysis/.
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
  if $PY "$script" > /tmp/repro_$(basename "$script" .py).log 2>&1; then
    echo "ran"
  else
    echo "FAILED TO RUN"
    sed 's/^/      /' /tmp/repro_$(basename "$script" .py).log | tail -12
    FAILED=1
  fi
}

# ---------------------------------------------------------------------------
# Stage 0 — restore the pre-net decomposition table.
# m1_decomposition.py emits arm-level values on the published (from-baseline)
# scale. v37_apply_net.py converts the Slentz arms to net-of-control. Running
# m1 first and v37 second reproduces the archived file exactly; the backup in
# analysis/ is the frozen input to v37 and doubles as the idempotency guard.
# ---------------------------------------------------------------------------
say "Stage 0 — baseline table"
run scripts/m1_decomposition.py

say "Stage 1 — net-scale recalibration"
run scripts/v37_apply_net.py

say "Stage 2 — meta-analysis (three axes)"
run scripts/m2_v313.py

say "Stage 3 — stratified by baseline glycaemic state"
run scripts/m2_stratified.py

say "Stage 4 — stratum x secretion interaction test"
run scripts/m4_interaction_test.py

say "Stage 5 — frozen arm-level dataset for figures"
run scripts/arms_v500.py

say "Stage 6 — weight-share traceability"
run scripts/m6_weight_shares.py

# ---------------------------------------------------------------------------
# Verification. Numeric outputs are compared byte-for-byte: the pipeline is
# deterministic, so any drift means the code and the archived result have
# diverged and the archived number is no longer reproducible.
# ---------------------------------------------------------------------------
say "Verification — regenerated outputs vs archived copies"

FILES=(
  analysis/decomposition.csv
  analysis/decomposition_v313.csv
  analysis/meta_A_v313.json
  analysis/stratified.json
  analysis/interaction_test_v314.json
  analysis/arms_v500.json
  analysis/weight_shares_v318.json
  analysis/v313_lee_secondphase.json
)

for f in "${FILES[@]}"; do
  printf '  %-42s' "$f"
  if [ ! -f "$f" ]; then
    echo "MISSING"
    FAILED=1
    continue
  fi
  if git -C "$ROOT" rev-parse --git-dir >/dev/null 2>&1; then
    if git -C "$ROOT" diff --quiet -- "$f"; then
      echo "PASS"
    else
      echo "FAIL (differs from committed copy)"
      FAILED=1
    fi
  else
    # No git history available: fall back to comparing against the manifest.
    want="$(grep -E "  $f\$" repro/expected_md5.txt | awk '{print $1}')"
    got="$(md5sum "$f" | awk '{print $1}')"
    if [ "$want" = "$got" ]; then echo "PASS"; else echo "FAIL (md5 $got, expected $want)"; FAILED=1; fi
  fi
done

# ---------------------------------------------------------------------------
# Spot-check the headline numbers so that a silent change in effect size
# cannot pass just because the files still hash-match each other.
# ---------------------------------------------------------------------------
say "Spot-check — headline effect sizes"

$PY - <<'PYEOF'
import json, sys
d = json.load(open("analysis/meta_A_v313.json"))["A_exercise_vs_control"]
expect = {
    "DI":  dict(k=8, pct=47.5, lo=16.9, hi=86.1),
    "SI":  dict(k=8, pct=42.7, lo=26.7, hi=60.7),
    "SEC": dict(k=7, pct=0.5,  lo=-17.6, hi=22.6),
}
bad = 0
for axis, e in expect.items():
    got = d[axis]
    ok = (got["k"] == e["k"]
          and round(got["pct"], 1) == e["pct"]
          and round(got["ci_lo_pct"], 1) == e["lo"]
          and round(got["ci_hi_pct"], 1) == e["hi"])
    print(f'  {axis:4s} k={got["k"]} {got["pct"]:+.1f}% '
          f'({got["ci_lo_pct"]:+.1f} to {got["ci_hi_pct"]:+.1f})  '
          f'{"PASS" if ok else "FAIL"}')
    bad += 0 if ok else 1
sys.exit(1 if bad else 0)
PYEOF
[ $? -ne 0 ] && FAILED=1

# ---------------------------------------------------------------------------
echo
if [ "$FAILED" -eq 0 ]; then
  echo "ALL CHECKS PASSED"
else
  echo "SOME CHECKS FAILED — see logs in /tmp/repro_*.log"
fi
exit "$FAILED"
