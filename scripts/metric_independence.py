#!/usr/bin/env python3
"""
Are the three "independent" prioritization metrics actually independent?

Reviewer point: ipTM and interface-pLDDT are both confidence outputs of the same
model and are probably correlated; Kd is a deterministic transform of dG. So
"all three metrics agree" may be closer to one or two effectively independent
signals. The manuscript already concedes the Kd/dG redundancy but not the
ipTM / i-pLDDT one.

This settles it from existing data. No new modelling.

Also reports the expected chance frequency of each ScanProsite motif under two
background models, as a first-pass check on whether 1,092 Tier 1 hits exceeds
what the motifs would return at random. This is NOT a substitute for a
composition-preserving shuffled-proteome control -- it is the back-of-envelope
that motivates running one.
"""
import os as _os
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import csv
from collections import defaultdict
import numpy as np
from scipy import stats

DATA = _os.path.join(_REPO, "data")
RES = _os.path.join(_REPO, "results")

rows = list(csv.DictReader(open(_os.path.join(DATA, "per_seed_metrics.csv"))))
METRICS = ["TCR_pMHC_ipTM", "interface_pLDDT", "dG_kcal_mol", "Kd"]

per = defaultdict(lambda: defaultdict(dict))
role = {}
for r in rows:
    if r["metric"] not in METRICS or r["provenance"] != "WORKBOOK":
        continue
    key = (r["system"], r["peptide_label"], r["sequence"])
    per[key][r["metric"]][int(r["seed"])] = float(r["value"])
    role[key] = r["role"]

cands = [k for k in per if role[k] == "candidate"]
print(f"Candidates: {len(cands)}\n")

def means(keys, m):
    return np.array([np.mean([per[k][m][s] for s in sorted(per[k][m])]) for k in keys])

# ---------------------------------------------------------------- 1. across candidates
print("=" * 72)
print("1. CORRELATION ACROSS CANDIDATE MEANS (n = 55)")
print("=" * 72)
labels = {"TCR_pMHC_ipTM": "ipTM", "interface_pLDDT": "i-pLDDT",
          "dG_kcal_mol": "dG", "Kd": "Kd"}
vals = {m: means(cands, m) for m in METRICS}
print(f"{'pair':26s} {'Pearson r':>10s} {'p':>10s} {'Spearman rho':>14s} {'p':>10s}")
for i, a in enumerate(METRICS):
    for b in METRICS[i + 1:]:
        r_p, p_p = stats.pearsonr(vals[a], vals[b])
        r_s, p_s = stats.spearmanr(vals[a], vals[b])
        print(f"{labels[a]+' vs '+labels[b]:26s} {r_p:10.3f} {p_p:10.2e} {r_s:14.3f} {p_s:10.2e}")

# ---------------------------------------------------------------- 2. within-candidate, across seeds
print()
print("=" * 72)
print("2. WITHIN-CANDIDATE SEED-LEVEL CORRELATION, ipTM vs i-pLDDT")
print("   (does a seed that scores high on one also score high on the other?)")
print("=" * 72)
rs = []
for k in cands:
    a = [per[k]["TCR_pMHC_ipTM"][s] for s in sorted(per[k]["TCR_pMHC_ipTM"])]
    b = [per[k]["interface_pLDDT"][s] for s in sorted(per[k]["interface_pLDDT"])]
    if len(set(a)) > 1 and len(set(b)) > 1:
        rs.append(stats.pearsonr(a, b)[0])
rs = np.array(rs)
print(f"candidates with variable values in both: {len(rs)} / {len(cands)}")
print(f"median within-candidate r = {np.median(rs):.3f}")
print(f"mean   within-candidate r = {np.mean(rs):.3f}")
print(f"fraction with r > 0       = {np.mean(rs > 0):.2%}")

# ---------------------------------------------------------------- 3. effective dimensionality
print()
print("=" * 72)
print("3. EFFECTIVE DIMENSIONALITY OF THE THREE PRIORITIZATION METRICS")
print("=" * 72)
three = ["TCR_pMHC_ipTM", "interface_pLDDT", "dG_kcal_mol"]
X = np.column_stack([vals[m] for m in three])
Z = (X - X.mean(0)) / X.std(0, ddof=1)
ev = np.linalg.eigvalsh(np.corrcoef(Z, rowvar=False))[::-1]
print("correlation-matrix eigenvalues:", np.round(ev, 3))
print(f"variance explained by PC1: {ev[0] / ev.sum():.1%}")
part = (ev.sum() ** 2) / np.sum(ev ** 2)
print(f"participation-ratio effective dimensionality: {part:.2f} of 3")

# Kd vs dG, to confirm the already-conceded redundancy on the same footing
r_kd = stats.pearsonr(vals["dG_kcal_mol"], np.log10(vals["Kd"]))[0]
print(f"\ndG vs log10(Kd) Pearson r = {r_kd:.4f}  (deterministic transform, as stated)")

# ---------------------------------------------------------------- 4. motif chance frequency
print()
print("=" * 72)
print("4. EXPECTED CHANCE FREQUENCY OF EACH ScanProsite MOTIF")
print("=" * 72)
# Swiss-Prot average amino-acid composition (fractional). Standard published
# background frequencies; used here only as a background approximation.
SP = dict(A=.0825, R=.0553, N=.0406, D=.0546, C=.0138, Q=.0393, E=.0672, G=.0707,
          H=.0227, I=.0591, L=.0965, K=.0580, M=.0241, F=.0386, P=.0470, S=.0656,
          T=.0534, W=.0110, Y=.0292, V=.0687)
MOTIFS = {
    "Tier1 #1": [["A","G","S","T"], ["L","I","V"], ["E","D"], ["A","G","S","T"], ["E","D"]],
    "Tier1 #2": [["A","G","S","T"], ["L","I","V"], ["A","G","S","T"], ["E","D"], ["L","I","V"]],
    "Tier1 #5": [["E","D"], ["L","I","V"], ["A","G","S","T"], ["L","I","V"], ["E","D"]],
    "Tier2 #3": [["L","I","V"], ["E","D"], ["A","G","S","T"], ["E","D"], ["E","D"], ["L","I","V"]],
    "Tier2 #4": [["L","I","V"], ["A","G","S","T"], ["L","I","V"], ["E","D"], ["L","I","V"]],
}
print(f"{'motif':12s} {'constrained':>12s} {'p(uniform)':>14s} {'1 in':>10s} "
      f"{'p(SwissProt)':>14s} {'1 in':>10s}")
for name, cls in MOTIFS.items():
    p_u = np.prod([len(c) / 20 for c in cls])
    p_s = np.prod([sum(SP[a] for a in c) for c in cls])
    print(f"{name:12s} {len(cls):12d} {p_u:14.3e} {1/p_u:10.0f} {p_s:14.3e} {1/p_s:10.0f}")

print("""
Interpretation: these are per-start-position probabilities under an i.i.d.
background. Multiplying by the number of scanned start positions in the searched
UniProtKB/Swiss-Prot subset gives the expected chance yield. That residue count
is not recorded in this archive, so the comparison against the observed 1,092
Tier 1 hits cannot be completed here -- which is precisely the gap. The proper
control is a composition-preserving shuffle of the same proteomes through the
identical ScanProsite pipeline, not this calculation.
""")

with open(_os.path.join(RES, "metric_independence.txt"), "w") as f:
    f.write("See console output of scripts/metric_independence.py\n")
print(f"(run: python3 scripts/metric_independence.py)")

# ---------------------------------------------------------------- 5. within-system
print()
print("=" * 72)
print("5. WITHIN-SYSTEM vs POOLED (pooling can create or hide correlation)")
print("=" * 72)
three = ["TCR_pMHC_ipTM", "interface_pLDDT", "dG_kcal_mol"]
print(f"{'group':26s} {'n':>4s} {'ipTM~ipLDDT':>13s} {'ipTM~dG':>10s} {'ipLDDT~dG':>11s} {'effdim':>7s}")
groups = [("POOLED (all systems)", cands)]
for s in ["HIP6/A2.11/DQ8", "HIP11/8.E3/DQ8-trans", "HIP11/8.E3/DQ8"]:
    groups.append((s, [k for k in cands if k[0] == s]))
for name, ks in groups:
    if len(ks) < 4:
        print(f"{name:26s} {len(ks):>4d}   (n too small)"); continue
    a, b, c = (means(ks, m) for m in three)
    Zg = np.column_stack([a, b, c]); Zg = (Zg - Zg.mean(0)) / Zg.std(0, ddof=1)
    e = np.linalg.eigvalsh(np.corrcoef(Zg, rowvar=False))[::-1]
    print(f"{name:26s} {len(ks):>4d} {stats.pearsonr(a,b)[0]:13.3f} "
          f"{stats.pearsonr(a,c)[0]:10.3f} {stats.pearsonr(b,c)[0]:11.3f} "
          f"{(e.sum()**2)/np.sum(e**2):7.2f}")
print("""
Pooling did BOTH things it can do: it HID the ipTM/i-pLDDT association (pooled
r = 0.560; within-system 0.648, 0.905, 0.972) and it CREATED the ipTM/dG one
(pooled r = -0.580; within-system +0.086, -0.240, -0.457). Within-system values
are the ones to report.
""")
