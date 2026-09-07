#!/usr/bin/env python3
"""
Phase 1.1 -- CSV spine with provenance tags.

Reads the ONLY located machine-readable source for the 55 Tier 1 candidates,
3 negative controls, and 3 TCRModel2-based native-HIP re-runs:
    ~/Downloads/Master Results Sheet.xlsx

Writes three tidy CSVs to Reproducibility_Pipeline/data/:
    per_seed_metrics.csv          -- one row per (system, peptide, seed, metric)
    per_peptide_summary.csv       -- one row per (system, peptide): RMSD, p-value,
                                      pass/fail, CI text, SD, as already stored in
                                      the workbook (not recomputed here)
    reference_complex_af3_metrics.csv -- the 5 reference systems' AlphaFold 3
                                      per-seed ipTM/pTM, from real job JSON files
                                      (genuinely RAW provenance)

Provenance levels (per the ground rules):
    RAW      = traceable to an actual model-output file on disk
    WORKBOOK = exists only as a pasted value in Master Results Sheet.xlsx;
               no upstream model-output artifact has been found
    MISSING  = not present anywhere; written explicitly as an empty row so the
               gap is visible in the data itself, not just in prose

No value in this script is invented, imputed, or estimated. Every WORKBOOK row
carries the exact workbook cell coordinates it came from.
"""
import csv
import json
import glob
import openpyxl
import os as _os
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))

XLSX = _os.path.join(_REPO, "data", "Master_Results_Sheet.xlsx")
DOWNLOADS = _os.path.join(_REPO, "alphafold3_reference_complexes")
OUT = _os.path.join(_REPO, "data")

wb = openpyxl.load_workbook(XLSX, data_only=True)
ws = wb.active
SHEET_NAME = "Master Results Sheet"

# ---------------------------------------------------------------------------
# System block row ranges, read directly from column A / B labels in the
# workbook (verified manually against the sheet before hardcoding here).
# ---------------------------------------------------------------------------
BLOCKS = [
    dict(system="HIP6/A2.11/DQ8",       header_row=1, cand_rows=range(3, 31),  native_row=31, control_row=32),
    dict(system="HIP11/8.E3/DQ8",       header_row=35, cand_rows=range(37, 44), native_row=44, control_row=45),
    dict(system="HIP11/8.E3/DQ8-trans", header_row=50, cand_rows=range(52, 72), native_row=72, control_row=73),
]

PER_SEED_METRIC_COLS = {
    "TCR_pMHC_ipTM":   list(range(4, 9)),    # cols D-H
    "interface_pLDDT": list(range(9, 14)),   # cols I-M
    "dG_kcal_mol":     list(range(14, 19)),  # cols N-R
    "Kd":              list(range(20, 25)),  # cols T-X
}
# columns known to NOT exist in this workbook at all, for any of the 61 peptides
MISSING_PER_SEED_METRICS = ["pLDDT", "pTM"]

SUMMARY_COLS = {
    "rmsd_angstrom": 26,       # Z
    "p_value_stored": 27,      # AA
    "pass_fail_stored": 28,    # AB
    "tcr_pmhc_iptm_ci_text": 29,  # AC
    "tcr_pmhc_iptm_sd_stored": 31,  # AE
    "model_selected_for_rmsd": 32,  # AF
}


def cell_ref(row, col):
    return ws.cell(row=row, column=col).coordinate


def get_peptide_info(row):
    label = ws.cell(row=row, column=1).value
    seq = ws.cell(row=row, column=2).value
    virus = ws.cell(row=row, column=3).value
    label = (label or "").strip()
    seq = (seq or "").strip() if seq else ""
    virus = (virus or "").strip() if virus else ""
    return label, seq, virus


def role_for(label, is_native, is_control):
    if is_native:
        return "native_HIP"
    if is_control:
        return "negative_control"
    return "candidate"


per_seed_rows = []
per_peptide_rows = []

for block in BLOCKS:
    system = block["system"]
    all_rows = list(block["cand_rows"]) + [block["native_row"], block["control_row"]]
    for row in all_rows:
        is_native = (row == block["native_row"])
        is_control = (row == block["control_row"])
        label, seq, virus = get_peptide_info(row)
        role = role_for(label, is_native, is_control)

        # --- per-seed metrics, present columns ---
        for metric, cols in PER_SEED_METRIC_COLS.items():
            for seed_idx, col in enumerate(cols, start=1):
                val = ws.cell(row=row, column=col).value
                provenance = "WORKBOOK"
                note = ""
                if val is None or (isinstance(val, str) and val.strip().upper() in ("N/A", "")):
                    provenance = "MISSING"
                    note = "cell present but blank/N/A in workbook"
                    val = ""
                per_seed_rows.append(dict(
                    system=system, peptide_label=label, sequence=seq, source_organism=virus,
                    role=role, seed=seed_idx, metric=metric, value=val,
                    provenance=provenance,
                    source_file=f"{XLSX}#{SHEET_NAME}!{cell_ref(row, col)}",
                    note=note,
                ))

        # --- per-seed metrics that do not exist anywhere in the workbook ---
        for metric in MISSING_PER_SEED_METRICS:
            for seed_idx in range(1, 6):
                per_seed_rows.append(dict(
                    system=system, peptide_label=label, sequence=seq, source_organism=virus,
                    role=role, seed=seed_idx, metric=metric, value="",
                    provenance="MISSING",
                    source_file="",
                    note="no pLDDT/pTM column exists in Master Results Sheet.xlsx for any of the 61 peptides",
                ))

        # --- per-peptide summary row ---
        summary = dict(system=system, peptide_label=label, sequence=seq, source_organism=virus, role=role)
        for key, col in SUMMARY_COLS.items():
            val = ws.cell(row=row, column=col).value
            summary[key] = val if val is not None else ""
            summary[f"{key}__cell"] = cell_ref(row, col)
        summary["provenance"] = "WORKBOOK"
        summary["source_file"] = f"{XLSX}#{SHEET_NAME}!row{row}"
        per_peptide_rows.append(summary)

# ---------------------------------------------------------------------------
# Write per_seed_metrics.csv
# ---------------------------------------------------------------------------
with open(f"{OUT}/per_seed_metrics.csv", "w", newline="") as f:
    fieldnames = ["system", "peptide_label", "sequence", "source_organism", "role",
                  "seed", "metric", "value", "provenance", "source_file", "note"]
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(per_seed_rows)

# ---------------------------------------------------------------------------
# Write per_peptide_summary.csv
# ---------------------------------------------------------------------------
with open(f"{OUT}/per_peptide_summary.csv", "w", newline="") as f:
    fieldnames = list(per_peptide_rows[0].keys())
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(per_peptide_rows)

# ---------------------------------------------------------------------------
# reference_complex_af3_metrics.csv -- genuinely RAW, from real AF3 job JSONs
# ---------------------------------------------------------------------------
AF3_JOBS = [
    ("HIP6/A2.11/DQ8",       "fold_hip6_tcr_a211_dq8"),
    ("HIP11/8.E3/DQ8",       "fold_hip_run_2"),
    ("HIP11/8.E3/DQ8-trans", "fold_hip11_gse8e3_dq8trans"),
    ("E2/HIP11/DQ2",         "fold_hip11_tcr_e2b_dq2"),
    ("A3.10/HIP6/DQ8",       "fold_hip6_tcr_a310_dq8"),
]

af3_rows = []
for system, folder in AF3_JOBS:
    for seed_idx in range(5):
        path = f"{DOWNLOADS}/{folder}/{folder}_summary_confidences_{seed_idx}.json"
        try:
            with open(path) as jf:
                d = json.load(jf)
            iptm = d.get("iptm")
            ptm = d.get("ptm")
            provenance = "RAW"
        except FileNotFoundError:
            iptm = ptm = ""
            provenance = "MISSING"
        af3_rows.append(dict(
            system=system, model_index=seed_idx, iptm=iptm, ptm=ptm,
            provenance=provenance, source_file=path,
        ))

with open(f"{OUT}/reference_complex_af3_metrics.csv", "w", newline="") as f:
    fieldnames = ["system", "model_index", "iptm", "ptm", "provenance", "source_file"]
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(af3_rows)

# ---------------------------------------------------------------------------
# Provenance summary report
# ---------------------------------------------------------------------------
from collections import Counter
c = Counter(r["provenance"] for r in per_seed_rows)
c_af3 = Counter(r["provenance"] for r in af3_rows)

print("=== per_seed_metrics.csv ===")
print(f"  total rows: {len(per_seed_rows)}")
for level in ["RAW", "WORKBOOK", "MISSING"]:
    print(f"  {level}: {c.get(level,0)}")

print()
print("=== per_peptide_summary.csv ===")
print(f"  total rows: {len(per_peptide_rows)}  (all WORKBOOK provenance)")

print()
print("=== reference_complex_af3_metrics.csv ===")
print(f"  total rows: {len(af3_rows)}")
for level in ["RAW", "WORKBOOK", "MISSING"]:
    print(f"  {level}: {c_af3.get(level,0)}")

print()
n_peptides = len(per_peptide_rows)
print(f"Peptides covered: {n_peptides} (should be 61: 28+1+1 HIP6, 7+1+1 cis, 20+1+1 trans)")
