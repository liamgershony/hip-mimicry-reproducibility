#!/usr/bin/env python3
"""
Statistical reanalysis using ONLY data already on disk. No new modelling.

Addresses the analysable half of Phase 4:
  - Mann-Whitney U instead of Welch's t (ipTM is bounded on [0,1], n=5, and
    normality is not established)
  - Benjamini-Hochberg FDR correction across all 55 candidate comparisons
  - Effect sizes (Cliff's delta) with percentile bootstrap CIs
  - An explicit statement of what n=5 can and cannot resolve

Reads: data/per_seed_metrics.csv
Writes: results/reanalysis_no_new_data.csv  +  console summary
"""
import csv, math, itertools
from collections import defaultdict
import numpy as np
from scipy import stats
import os as _os
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))

DATA = _os.path.join(_REPO, "data")
RES = _os.path.join(_REPO, "results")

rows = list(csv.DictReader(open(f"{DATA}/per_seed_metrics.csv")))
g = defaultdict(dict)
meta = {}
for r in rows:
    if r["metric"] != "TCR_pMHC_ipTM" or r["provenance"] != "WORKBOOK":
        continue
    key = (r["system"], r["peptide_label"], r["sequence"])
    g[key][int(r["seed"])] = float(r["value"])
    meta[key] = r["role"]

controls, candidates = {}, {}
for key, seeds in g.items():
    vals = [seeds[s] for s in sorted(seeds)]
    if meta[key] == "negative_control":
        controls[key[0]] = vals
    elif meta[key] == "candidate":
        candidates[key] = vals

assert len(candidates) == 55 and len(controls) == 3

def cliffs_delta(a, b):
    """P(a>b) - P(a<b); +1 = complete separation above, -1 = below."""
    gt = sum(1 for x, y in itertools.product(a, b) if x > y)
    lt = sum(1 for x, y in itertools.product(a, b) if x < y)
    return (gt - lt) / (len(a) * len(b))

def boot_ci_delta(a, b, n=10000, seed=0):
    rng = np.random.default_rng(seed)
    a, b = np.array(a), np.array(b)
    ds = [cliffs_delta(rng.choice(a, len(a), replace=True),
                        rng.choice(b, len(b), replace=True)) for _ in range(n)]
    return float(np.percentile(ds, 2.5)), float(np.percentile(ds, 97.5))

# ---- minimum achievable p for a two-sided MWU with n1=n2=5 ----
best = stats.mannwhitneyu([10,11,12,13,14], [1,2,3,4,5], alternative="two-sided").pvalue
print(f"Floor on two-sided Mann-Whitney p with n=5 vs n=5 (perfect separation): p = {best:.4f}")
print(f"  -> no candidate can score below this, however large the true effect.\n")

out = []
for key, vals in candidates.items():
    system, label, seq = key
    ctrl = controls[system]
    w_p = stats.ttest_ind(vals, ctrl, equal_var=False).pvalue
    u_stat, u_p = stats.mannwhitneyu(vals, ctrl, alternative="two-sided")
    d = cliffs_delta(vals, ctrl)
    lo, hi = boot_ci_delta(vals, ctrl)
    out.append(dict(system=system, peptide=label.strip(), sequence=seq,
                     mean_iptm=round(float(np.mean(vals)), 4),
                     control_mean=round(float(np.mean(ctrl)), 4),
                     welch_p=round(float(w_p), 6),
                     mannwhitney_p=round(float(u_p), 6),
                     cliffs_delta=round(d, 3),
                     delta_ci_low=round(lo, 3), delta_ci_high=round(hi, 3)))

# ---- Benjamini-Hochberg across all 55 ----
def bh(pvals, alpha=0.05):
    m = len(pvals)
    order = np.argsort(pvals)
    ranked = np.array(pvals)[order]
    thresh = (np.arange(1, m + 1) / m) * alpha
    passed = ranked <= thresh
    k = np.max(np.where(passed)[0]) + 1 if passed.any() else 0
    reject = np.zeros(m, dtype=bool)
    if k:
        reject[order[:k]] = True
    # adjusted p-values (step-up, monotone)
    adj = np.minimum.accumulate((ranked * m / np.arange(1, m + 1))[::-1])[::-1]
    adj_full = np.empty(m)
    adj_full[order] = np.minimum(adj, 1.0)
    return reject, adj_full

for label, key in [("welch", "welch_p"), ("mannwhitney", "mannwhitney_p")]:
    ps = [r[key] for r in out]
    rej, adj = bh(ps)
    for r, a, j in zip(out, adj, rej):
        r[f"{label}_p_BH"] = round(float(a), 6)
        r[f"{label}_sig_after_BH"] = bool(j)

with open(f"{RES}/reanalysis_no_new_data.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
    w.writeheader(); w.writerows(out)

n_w_raw = sum(1 for r in out if r["welch_p"] < 0.05)
n_w_bh = sum(1 for r in out if r["welch_sig_after_BH"])
n_u_raw = sum(1 for r in out if r["mannwhitney_p"] < 0.05)
n_u_bh = sum(1 for r in out if r["mannwhitney_sig_after_BH"])

print("Significant out of 55 candidates:")
print(f"  Welch's t, uncorrected         : {n_w_raw}   <- the published '48/55'")
print(f"  Welch's t, BH-corrected        : {n_w_bh}")
print(f"  Mann-Whitney U, uncorrected    : {n_u_raw}")
print(f"  Mann-Whitney U, BH-corrected   : {n_u_bh}")
print()

dsep = [r for r in out if r["cliffs_delta"] == 1.0]
print(f"Candidates with complete rank separation from their control "
      f"(Cliff's delta = +1.0): {len(dsep)} / 55")
neg = [r for r in out if r["cliffs_delta"] < 0]
print(f"Candidates scoring BELOW their control (delta < 0): {len(neg)}")
if neg:
    for r in neg:
        print(f"    {r['sequence']:16s} {r['system']:22s} delta={r['cliffs_delta']:+.2f} "
              f"mean={r['mean_iptm']} vs ctrl {r['control_mean']}")
print()

lead = [r for r in out if r["sequence"] == "RRNVATLQAENVTG"][0]
print("Lead candidate RRNVATLQAENVTG vs its own negative control:")
for k in ["mean_iptm", "control_mean", "welch_p", "welch_p_BH", "mannwhitney_p",
          "mannwhitney_p_BH", "cliffs_delta", "delta_ci_low", "delta_ci_high"]:
    print(f"    {k:20s} {lead[k]}")
print()
print(f"Wrote {RES}/reanalysis_no_new_data.csv")
