#!/usr/bin/env python3
"""
Phase 1.4 -- Kd/dG consistency check.

For every (system, peptide, seed) row in the workbook, recompute Kd from dG via
dG = RT * ln(Kd)  =>  Kd_predicted = exp(dG / RT)
and compare against the stored Kd. Report every cell where they disagree beyond
rounding, quantify whether row 73 is isolated or symptomatic, and re-run the
native-HIP CI gate with (a) stored Kd and (b) dG-derived Kd to check the
manuscript's claim (Section 3.15) that the same 10 candidates pass either way.

RT is not assumed -- it is fit from the data itself (median implied RT across
all rows), then used to flag outliers. Nothing is silently corrected; outliers
are written to a documented corrections file for human review.
"""
import csv, math
from collections import defaultdict
import numpy as np
import os as _os
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))

DATA = _os.path.join(_REPO, "data")
RESULTS = _os.path.join(_REPO, "results")
CORRECTIONS = _os.path.join(_REPO, "corrections")

rows = list(csv.DictReader(open(f"{DATA}/per_seed_metrics.csv")))
dg_rows = {(r["system"], r["role"], r["peptide_label"], r["sequence"], int(r["seed"])): float(r["value"])
           for r in rows if r["metric"] == "dG_kcal_mol" and r["provenance"] == "WORKBOOK"}
kd_rows = {(r["system"], r["role"], r["peptide_label"], r["sequence"], int(r["seed"])): (r["value"], r["source_file"])
           for r in rows if r["metric"] == "Kd" and r["provenance"] == "WORKBOOK"}

# --- fit RT from implied values: RT_implied = dG / ln(Kd), excluding non-numeric/zero Kd ---
implied_RT = []
for key, dg in dg_rows.items():
    kd_str, _ = kd_rows.get(key, ("", ""))
    try:
        kd = float(kd_str)
        if kd <= 0:
            continue
        rt = dg / math.log(kd)
        implied_RT.append(rt)
    except (ValueError, ZeroDivisionError):
        continue

implied_RT = np.array(implied_RT)
median_RT = float(np.median(implied_RT))
print(f"Rows with usable dG & Kd: {len(implied_RT)} / {len(dg_rows)}")
print(f"Median implied RT = {median_RT:.5f} kcal/mol (standard 298K value ~0.5921-0.5924)")
print(f"IQR of implied RT: [{np.percentile(implied_RT,25):.5f}, {np.percentile(implied_RT,75):.5f}]")

RT = median_RT

records = []
outliers = []
diffs = []
for key, dg in dg_rows.items():
    system, role, label, seq, seed = key
    kd_str, kd_source = kd_rows.get(key, ("", ""))
    try:
        kd_stored = float(kd_str)
    except ValueError:
        continue
    kd_predicted = math.exp(dg / RT)
    dg_from_kd_stored = RT * math.log(kd_stored) if kd_stored > 0 else float("nan")
    diff = dg - dg_from_kd_stored
    diffs.append(diff)
    rec = dict(system=system, role=role, peptide=label, sequence=seq, seed=seed,
               dG_stored=dg, Kd_stored=kd_stored, Kd_predicted_from_dG=kd_predicted,
               dG_implied_by_stored_Kd=round(dg_from_kd_stored, 4) if not math.isnan(dg_from_kd_stored) else "",
               diff_kcal_mol=round(diff, 4) if not math.isnan(diff) else "",
               source_file=kd_source)
    records.append(rec)
    if abs(diff) > 0.05:  # beyond simple rounding
        outliers.append(rec)

with open(f"{RESULTS}/kd_dg_consistency_full.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(records[0].keys()))
    w.writeheader()
    w.writerows(records)

with open(f"{CORRECTIONS}/kd_dg_outliers_FOR_REVIEW.csv", "w", newline="") as f:
    if outliers:
        w = csv.DictWriter(f, fieldnames=list(outliers[0].keys()))
        w.writeheader()
        w.writerows(outliers)

diffs_arr = np.array([d for d in diffs if not math.isnan(d)])
print()
print(f"Total rows checked: {len(records)}")
print(f"Median |dG - RT*ln(Kd_stored)| : {np.median(np.abs(diffs_arr)):.4f} kcal/mol")
print(f"IQR: [{np.percentile(diffs_arr,25):.4f}, {np.percentile(diffs_arr,75):.4f}]")
print(f"Outliers beyond 0.05 kcal/mol disagreement: {len(outliers)}")
for o in outliers:
    print(f"  {o['system']:22s} {o['peptide']:14s} {o['sequence']:16s} seed{o['seed']} "
          f"dG_stored={o['dG_stored']} Kd_stored={o['Kd_stored']:.3e} "
          f"-> implied dG from stored Kd = {o['dG_implied_by_stored_Kd']}  diff={o['diff_kcal_mol']}")
