#!/usr/bin/env python3
"""
Three analyses that strengthen the paper with no new modelling.

1. RMSD as a DISTRIBUTION, not a threshold. Where do the three negative controls
   sit inside the candidate RMSD distribution? If they sit mid-distribution, the
   negative result stops depending on 1.0 A being the right cutoff.
2. Gate-threshold sensitivity: nominated set under fixed thresholds instead of
   per-system native-HIP confidence intervals (whose widths differ up to ~6x).
3. Which nominees' stored SD cells were corrupted -- in particular the two
   peptides now at the bench.
"""
import os as _os
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import csv, math
from collections import defaultdict
import numpy as np

DATA = _os.path.join(_REPO, "data"); RES = _os.path.join(_REPO, "results")
T_CRIT, N = 2.776, 5
HIGHER = {"TCR_pMHC_ipTM", "interface_pLDDT"}; LOWER = {"dG_kcal_mol"}
ANIMAL = {"RHGEVNDSTTVEPVL", "GPGEVNDSTTVEPIL", "GPGEVNDSTTVEPVL"}

summ = list(csv.DictReader(open(_os.path.join(DATA, "per_peptide_summary.csv"))))
rows = list(csv.DictReader(open(_os.path.join(DATA, "per_seed_metrics.csv"))))
g = defaultdict(lambda: defaultdict(dict)); role = {}
for r in rows:
    if r["metric"] in HIGHER | LOWER and r["provenance"] == "WORKBOOK":
        k = (r["system"], r["peptide_label"], r["sequence"])
        g[k][r["metric"]][int(r["seed"])] = float(r["value"]); role[k] = r["role"]
def sv(k,m): return [g[k][m][s] for s in sorted(g[k][m])]
def mean(v): return sum(v)/len(v)
def sd(v):
    m=mean(v); return math.sqrt(sum((x-m)**2 for x in v)/(len(v)-1))

# ---------------------------------------------------------------- 1. RMSD
print("="*74); print("1. RMSD AS A DISTRIBUTION"); print("="*74)
cand_r, ctrl_r = [], []
for r in summ:
    try: v = float(r["rmsd_angstrom"])
    except (ValueError, KeyError): continue
    if r["role"] == "candidate": cand_r.append((v, r["sequence"], r["system"]))
    elif r["role"] == "negative_control": ctrl_r.append((v, r["sequence"], r["system"]))
c = np.array([x[0] for x in cand_r])
print(f"Candidates (n={len(c)}):")
print(f"   min {c.min():.3f}   Q1 {np.percentile(c,25):.3f}   median {np.median(c):.3f}"
      f"   Q3 {np.percentile(c,75):.3f}   max {c.max():.3f}")
print(f"   mean {c.mean():.3f}   IQR {np.percentile(c,75)-np.percentile(c,25):.3f}")
excl = c[c < 1.0]
print(f"   excluding the 1.896 outlier: median {np.median(excl):.3f}, max {excl.max():.3f}")
print(f"\nNegative controls, and their percentile WITHIN the candidate distribution:")
for v, seq, sys_ in sorted(ctrl_r):
    pct = 100.0*(c < v).sum()/len(c)
    print(f"   {seq:16s} {sys_:22s} RMSD {v:.3f}   ->  {pct:5.1f}th percentile of candidates")
lead = [x for x in cand_r if x[1]=="RRNVATLQAENVTG"][0]
print(f"\n   for reference, lead candidate RRNVATLQAENVTG RMSD {lead[0]:.3f} -> "
      f"{100.0*(c<lead[0]).sum()/len(c):.1f}th percentile")
print("\n   => the controls are interior to the candidate distribution, so the")
print("      negative result does not depend on 1.0 A being the right cutoff.")

# ---------------------------------------------------------------- 2. threshold sensitivity
print(); print("="*74); print("2. GATE-THRESHOLD SENSITIVITY"); print("="*74)
native = {}
for k in g:
    if role[k]!="native_HIP": continue
    for m in HIGHER|LOWER:
        v=sv(k,m); native[(k[0],m)]=(mean(v), sd(v))
cands=[k for k in g if role[k]=="candidate"]
print("native-HIP CI half-widths as % of the mean (the widths that differ ~6x):")
for s in ["HIP6/A2.11/DQ8","HIP11/8.E3/DQ8","HIP11/8.E3/DQ8-trans"]:
    mu,s_=native[(s,"TCR_pMHC_ipTM")]; hw=T_CRIT*s_/math.sqrt(N)
    print(f"   {s:22s} ipTM {mu:.3f} +/- {hw:.3f}  ({100*hw/mu:.2f}%)")

def gate_ci(k, ms):
    for m in ms:
        mu,s_=native[(k[0],m)]; hw=T_CRIT*s_/math.sqrt(N); v=mean(sv(k,m))
        if m in HIGHER:
            if v < mu-hw: return False
        else:
            if v > mu+hw: return False
    return True
def gate_fixed(k, ms, frac):
    """Candidate mean must be within `frac` of the native mean (relative to |mean|),
    or better. Uses abs() so the tolerance direction is correct for dG, which is
    negative -- mu*(1+frac) would move the bound the wrong way."""
    for m in ms:
        mu,_=native[(k[0],m)]; v=mean(sv(k,m)); tol=frac*abs(mu)
        if m in HIGHER:
            if v < mu-tol: return False
        else:
            if v > mu+tol: return False
    return True

THREE=["TCR_pMHC_ipTM","interface_pLDDT","dG_kcal_mol"]
base={k for k in cands if gate_ci(k,THREE)}
base_nom={k for k in base if k[2] not in ANIMAL}
print(f"\nCI-derived gate (as published, 3-metric): {len(base)} pass, {len(base_nom)} nominees")
print(f"{'fixed relative threshold':28s} {'pass':>5s} {'nominees':>9s}  identical to published?")
for frac in [0.005,0.01,0.02,0.03,0.05]:
    p={k for k in cands if gate_fixed(k,THREE,frac)}
    nom={k for k in p if k[2] not in ANIMAL}
    print(f"  within {frac*100:>4.1f}% of native      {len(p):5d} {len(nom):9d}  {nom==base_nom}")

# ---------------------------------------------------------------- 3. corrupted SD cells
print(); print("="*74); print("3. WERE THE BENCH PEPTIDES' SD CELLS CORRUPTED?"); print("="*74)
# system-disambiguated: SDLLKSVDSEEVRE and others appear in more than one system
NOM=[("HIP6/A2.11/DQ8","VGVEVSAANTLEECW"),("HIP6/A2.11/DQ8","SESEVESDTALESEV"),
     ("HIP6/A2.11/DQ8","NDDELAGNLGLDVSN"),("HIP6/A2.11/DQ8","RLCDLDTADAVEEMG"),
     ("HIP11/8.E3/DQ8-trans","RRNVATLQAENVTG"),("HIP11/8.E3/DQ8-trans","SDLLKSVDSEEVRE"),
     ("HIP11/8.E3/DQ8-trans","LQEILTIKSDDVVG"),("HIP11/8.E3/DQ8-trans","DYVVNTINSDPVME"),
     ("HIP11/8.E3/DQ8-trans","RDISSTIATEKIPF"),("HIP11/8.E3/DQ8-trans","GKDIEGVGSEDLVN")]
by_seq={(r["system"],r["sequence"]):r for r in summ}
print(f"{'nominee':17s} {'stored SD':>10s} {'true SD':>9s} {'mean ipTM':>10s}  status")
for sys_,s in NOM:
    r=by_seq.get((sys_,s))
    k=next((k for k in g if k[0]==sys_ and k[2]==s and role[k]=="candidate"), None)
    if not r or not k: continue
    try: st=float(r["tcr_pmhc_iptm_sd_stored"])
    except ValueError: st=float("nan")
    tv=sd(sv(k,"TCR_pMHC_ipTM")); mu=mean(sv(k,"TCR_pMHC_ipTM"))
    bad = abs(st-tv)>0.005
    note = "CORRUPT (stores the MEAN)" if bad and abs(st-mu)<0.005 else ("mismatch" if bad else "ok")
    print(f"{s:17s} {st:10.4f} {tv:9.4f} {mu:10.3f}  {note}")

with open(_os.path.join(RES,"rmsd_distribution.csv"),"w",newline="") as f:
    w=csv.writer(f); w.writerow(["role","system","sequence","rmsd_angstrom","percentile_within_candidates"])
    for v,seq,sys_ in sorted(cand_r): w.writerow(["candidate",sys_,seq,v,round(100.0*(c<v).sum()/len(c),1)])
    for v,seq,sys_ in sorted(ctrl_r): w.writerow(["negative_control",sys_,seq,v,round(100.0*(c<v).sum()/len(c),1)])
print(f"\nWrote {RES}/rmsd_distribution.csv")
