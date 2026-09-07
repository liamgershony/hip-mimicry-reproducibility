#!/bin/bash
# One-command regeneration of every derived table and statistic in this archive.
# Requires: python3 with openpyxl, scipy, numpy.
# (rebuild_figure3.py additionally needs `tectonic` to compile the TikZ figure;
#  it is skipped automatically if tectonic is not installed.)
set -e
cd "$(dirname "$0")"

echo "1/7  Building CSV spine from data/Master_Results_Sheet.xlsx ..."
python3 scripts/build_csv_spine.py
echo
echo "2/7  Recomputing Welch's t-test (the 48/55 claim) ..."
python3 scripts/recompute_welch.py
echo
echo "3/7  Recomputing native-HIP CI gate (12 pass -> 10 nominated) ..."
python3 scripts/recompute_ci_gate.py
echo
echo "4/7  Kd/dG internal consistency across all 305 seed values ..."
python3 scripts/kd_dg_consistency.py
echo
echo "5/7  Sensitivity reanalysis: Mann-Whitney + Benjamini-Hochberg + Cliff's delta ..."
python3 scripts/reanalysis_no_new_data.py
echo
echo "6/7  Deriving the pre-declared reproduction tolerance band ..."
python3 scripts/derive_tolerance_band.py
echo
if command -v tectonic >/dev/null 2>&1; then
  echo "7/7  Rebuilding Figure 3 from verified data ..."
  python3 scripts/rebuild_figure3.py
else
  echo "7/7  SKIPPED (tectonic not installed) - Figure 3 rebuild"
fi
echo
echo "Done. Outputs in results/ ; flagged data issues in corrections/."
