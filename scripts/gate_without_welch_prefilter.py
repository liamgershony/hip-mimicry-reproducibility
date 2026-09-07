#!/usr/bin/env python3
"""
Does the nominated set depend on the confounded Welch step?

The manuscript applies the four-metric native-HIP CI gate only to the 48
candidates that passed the Welch-vs-control comparison. But that comparison is
now disclosed as selection-confounded. If the nominated set descends from a
filter the paper disowns, that is a real structural problem.

This applies the identical four-metric gate to ALL 55 candidates and compares.
No new modelling; runs entirely on the existing records.
"""
import os as _os
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import csv, math
from collections import defaultdict

DATA = _os.path.join(_REPO, "data")
RES = _os.path.join(_REPO, "results")
T_CRIT, N = 2.776, 5

HIGHER_BETTER = {"TCR_pMHC_ipTM", "interface_pLDDT"}
LOWER_BETTER = {"dG_kcal_mol", "Kd"}
METRICS = HIGHER_BETTER | LOWER_BETTER

# animal-associated strains excluded by the human-association nomination filter
ANIMAL = {"RHGEVNDSTTVEPVL", "GPGEVNDSTTVEPIL", "GPGEVNDSTTVEPVL"}

rows = list(csv.DictReader(open(_os.path.join(DATA, "per_seed_metrics.csv"))))
g = defaultdict(lambda: defaultdict(dict))
role = {}
for r in rows:
    if r["metric"] not in METRICS or r["provenance"] != "WORKBOOK":
        continue
    k = (r["system"], r["peptide_label"], r["sequence"])
    g[k][r["metric"]][int(r["seed"])] = float(r["value"])
    role[k] = r["role"]

def mean(v): return sum(v) / len(v)
def sd(v):
    m = mean(v); return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))
def seeds(k, m): return [g[k][m][s] for s in sorted(g[k][m])]

# native-HIP CI per (system, metric)
native_ci = {}
for k in g:
    if role[k] != "native_HIP": continue
    for m in METRICS:
        v = seeds(k, m); mu = mean(v); hw = T_CRIT * sd(v) / math.sqrt(N)
        native_ci[(k[0], m)] = (mu - hw, mu + hw)

def passes_gate(k):
    for m in METRICS:
        mu = mean(seeds(k, m)); lo, hi = native_ci[(k[0], m)]
        if m in HIGHER_BETTER:
            if mu < lo: return False
        else:
            if mu > hi: return False
    return True

cands = [k for k in g if role[k] == "candidate"]
assert len(cands) == 55, len(cands)

welch = {(r["system"], r["peptide"], r["sequence"])
         for r in csv.DictReader(open(_os.path.join(RES, "welch_recompute.csv")))
         if r["recomputed_pass"] == "True"}

all55 = {k for k in cands if passes_gate(k)}
from48 = {k for k in cands if k in welch and passes_gate(k)}

print(f"Candidates: {len(cands)}   Welch-significant: {len(welch)}\n")
print(f"Four-metric gate applied to the 48 Welch-significant : {len(from48)} pass")
print(f"Four-metric gate applied to ALL 55 candidates         : {len(all55)} pass")
print(f"IDENTICAL SET: {all55 == from48}\n")

extra = all55 - from48
if extra:
    print("Candidates that pass the gate but were excluded by the Welch step:")
    for k in sorted(extra):
        print(f"  {k[2]:16s} {k[0]:22s} mean ipTM={mean(seeds(k,'TCR_pMHC_ipTM')):.3f}")
else:
    print("No candidate excluded by the Welch step would have passed the gate.")

# what the seven non-significant candidates actually look like
print("\nThe 7 candidates the Welch step excluded:")
for k in sorted(k for k in cands if k not in welch):
    mu = mean(seeds(k, "TCR_pMHC_ipTM"))
    lo, _ = native_ci[(k[0], "TCR_pMHC_ipTM")]
    print(f"  {k[2]:16s} {k[0]:22s} mean ipTM={mu:.3f}  "
          f"(native CI lower bound {lo:.3f}) -> gate {'PASS' if passes_gate(k) else 'FAIL'}")

nom_all = {k for k in all55 if k[2] not in ANIMAL}
nom_48 = {k for k in from48 if k[2] not in ANIMAL}
print(f"\nAfter the human-association filter:")
print(f"  from the 48 : {len(nom_48)} nominees")
print(f"  from all 55 : {len(nom_all)} nominees")
print(f"  IDENTICAL: {nom_all == nom_48}")

with open(_os.path.join(RES, "gate_without_welch_prefilter.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["system", "peptide", "sequence", "mean_iptm",
                "welch_significant", "passes_four_metric_gate", "animal_associated"])
    for k in sorted(cands):
        w.writerow([k[0], k[1].strip(), k[2], round(mean(seeds(k, "TCR_pMHC_ipTM")), 4),
                    k in welch, passes_gate(k), k[2] in ANIMAL])
print(f"\nWrote {RES}/gate_without_welch_prefilter.csv")
