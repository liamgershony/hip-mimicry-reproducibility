#!/usr/bin/env python3
"""
Phase 1.2a -- recompute the 48/55 Welch's t-test result from per_seed_metrics.csv
(sourced from the actual workbook, not a screenshot transcription).

Two-sided unpaired Welch's t-test, TCR-pMHC ipTM, candidate's 5 seeds vs. its
own system's negative control's 5 seeds. alpha = 0.05. Matches manuscript Eq. 1.
"""
import csv
from collections import defaultdict
from scipy import stats
import os as _os
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))

DATA = _os.path.join(_REPO, "data")

rows = list(csv.DictReader(open(f"{DATA}/per_seed_metrics.csv")))
iptm_rows = [r for r in rows if r["metric"] == "TCR_pMHC_ipTM"]

# group into {(system, peptide_label): {seed: value}}
grouped = defaultdict(dict)
for r in iptm_rows:
    key = (r["system"], r["peptide_label"], r["role"], r["sequence"])
    grouped[key][int(r["seed"])] = float(r["value"])

controls = {}
candidates = {}
for (system, label, role, seq), seeds in grouped.items():
    vals = [seeds[s] for s in sorted(seeds)]
    if role == "negative_control":
        controls[system] = vals
    elif role == "candidate":
        candidates[(system, label, seq)] = vals

assert len(candidates) == 55, f"expected 55 candidates, got {len(candidates)}"
assert len(controls) == 3, f"expected 3 controls, got {len(controls)}"

summary = list(csv.DictReader(open(f"{DATA}/per_peptide_summary.csv")))
stored_p = {(r["system"], r["peptide_label"], r["sequence"]): r["p_value_stored"] for r in summary}
stored_pass = {(r["system"], r["peptide_label"], r["sequence"]): r["pass_fail_stored"] for r in summary}

results = []
n_sig = 0
flips = []
for (system, label, seq), vals in candidates.items():
    ctrl = controls[system]
    t, p = stats.ttest_ind(vals, ctrl, equal_var=False)
    new_pass = p < 0.05
    if new_pass:
        n_sig += 1
    sp_raw = stored_p.get((system, label, seq), "")
    sp_pass = stored_pass.get((system, label, seq), "")
    try:
        sp = float(sp_raw)
        stored_sig = sp < 0.05
    except ValueError:
        sp = None
        stored_sig = (sp_pass == "PASS")
    flip = (new_pass != stored_sig)
    results.append(dict(system=system, peptide=label, sequence=seq,
                         recomputed_p=round(p, 6), recomputed_pass=new_pass,
                         stored_p=sp_raw, stored_pass=sp_pass, flip=flip))
    if flip:
        flips.append((system, label, seq, sp_raw, sp_pass, p, new_pass))

with open(_os.path.join(_REPO, "results", "welch_recompute.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
    w.writeheader()
    w.writerows(results)

print(f"Total candidates: {len(results)}")
print(f"Recomputed significant (p<0.05): {n_sig} / 55")
print(f"Manuscript claims: 48 / 55")
print(f"MATCH: {n_sig == 48}")
print()
print("Flips vs. workbook's own stored Pass/Fail column:")
for f_ in flips:
    print(f"  {f_}")
print(f"Total flips: {len(flips)}")
