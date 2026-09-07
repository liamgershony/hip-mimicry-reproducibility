#!/usr/bin/env python3
"""
Phase 1.2b -- recompute the native-HIP 95% CI gate (Eq. 2) and confirm the
10-candidate / 4-0-6 split claimed in the manuscript.

Candidate means must fall within, or perform better than, the native-HIP 95%
CI on ALL FOUR metrics: TCR-pMHC ipTM (higher better), interface-pLDDT (higher
better), dG (lower/more negative better), Kd (lower better).

t_crit = 2.776 (two-sided 95%, df=4, n=5), per manuscript Eq. 2.

This only evaluates candidates that already passed the Welch's t-test gate
(48/55, per recompute_welch.py), matching the manuscript's stated pipeline
order (Welch gate first, native-HIP CI gate second, §3.6-3.7).
"""
import csv, math
from collections import defaultdict
import os as _os
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))

DATA = _os.path.join(_REPO, "data")
RESULTS = _os.path.join(_REPO, "results")

T_CRIT = 2.776
N = 5

rows = list(csv.DictReader(open(f"{DATA}/per_seed_metrics.csv")))
METRICS_HIGHER_BETTER = {"TCR_pMHC_ipTM", "interface_pLDDT"}
METRICS_LOWER_BETTER = {"dG_kcal_mol", "Kd"}
METRICS = METRICS_HIGHER_BETTER | METRICS_LOWER_BETTER

# group: {(system, role, label, seq, metric): [5 values]}
grouped = defaultdict(dict)
for r in rows:
    if r["metric"] not in METRICS or r["provenance"] != "WORKBOOK":
        continue
    key = (r["system"], r["role"], r["peptide_label"], r["sequence"], r["metric"])
    grouped[key][int(r["seed"])] = float(r["value"])

def mean_sd(vals):
    n = len(vals)
    m = sum(vals) / n
    var = sum((v - m) ** 2 for v in vals) / (n - 1)
    return m, math.sqrt(var)

def ci(vals):
    m, sd = mean_sd(vals)
    hw = T_CRIT * sd / math.sqrt(N)
    return m, m - hw, m + hw

# native CI per system per metric
native_ci = {}
for (system, role, label, seq, metric), seeds in grouped.items():
    if role != "native_HIP":
        continue
    vals = [seeds[s] for s in sorted(seeds)]
    m, lo, hi = ci(vals)
    native_ci[(system, metric)] = dict(mean=m, lo=lo, hi=hi)

# load Welch-significant candidates (must match recompute_welch.py output)
welch = list(csv.DictReader(open(f"{RESULTS}/welch_recompute.csv")))
sig_candidates = {(r["system"], r["peptide"], r["sequence"]) for r in welch if r["recomputed_pass"] == "True"}
print(f"Welch-significant candidates carried into CI gate: {len(sig_candidates)}")

# evaluate each significant candidate on all 4 metrics
report = []
passing = []
for (system, role, label, seq, metric), seeds in grouped.items():
    pass  # placeholder, real loop below

candidate_metrics = defaultdict(dict)
for (system, role, label, seq, metric), seeds in grouped.items():
    if role != "candidate":
        continue
    if (system, label, seq) not in sig_candidates:
        continue
    vals = [seeds[s] for s in sorted(seeds)]
    m = sum(vals) / len(vals)
    candidate_metrics[(system, label, seq)][metric] = m

for (system, label, seq), metric_means in candidate_metrics.items():
    row = dict(system=system, peptide=label, sequence=seq)
    all_pass = True
    for metric in METRICS:
        m = metric_means.get(metric)
        nci = native_ci.get((system, metric))
        if m is None or nci is None:
            row[f"{metric}_result"] = "MISSING"
            all_pass = False
            continue
        if metric in METRICS_HIGHER_BETTER:
            ok = m >= nci["lo"]
        else:
            ok = m <= nci["hi"]
        row[f"{metric}_mean"] = round(m, 4)
        row[f"{metric}_native_ci"] = f"[{nci['lo']:.4f}, {nci['hi']:.4f}]"
        row[f"{metric}_result"] = "PASS" if ok else "FAIL"
        if not ok:
            all_pass = False
    row["overall_pass_all_4"] = all_pass
    report.append(row)
    if all_pass:
        passing.append((system, label, seq))

with open(f"{RESULTS}/ci_gate_recompute.csv", "w", newline="") as f:
    fieldnames = list(report[0].keys())
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(report)

print(f"\nCandidates passing all 4 metrics (recomputed, using raw Kd for the Kd column): {len(passing)}")
by_system = defaultdict(int)
for system, label, seq in passing:
    by_system[system] += 1
    print(f"  {system:22s} {label:14s} {seq}")
print()
print("Split by system:", dict(by_system))
print("Manuscript claims: 10 total, split 4 / 0 / 6")
print("MATCH:", len(passing) == 10)
