#!/usr/bin/env python3
"""
Two things that improve rigour with no new modelling.

PART 1 -- The four-metric gate is not four-dimensional.
  Kd = exp(dG/RT) exactly, so Kd and dG are one measurement. ipTM and
  interface-pLDDT correlate at 0.905 in the trans system. The paper discloses
  this but the gate was never redesigned. Recompute the nominated set under
  4-metric, 3-metric (drop Kd), and 2-metric (ipTM + dG) definitions and report
  whether the answer changes.

PART 2 -- Systematic audit of every stored workbook value against recomputation.
  A single stored/recomputed p-value discrepancy is already known
  (0.0658 -> 0.0404). The open question is whether there are others nobody
  checked. This audits every stored column that can be independently derived.
"""
import os as _os
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import csv, math
from collections import defaultdict
import numpy as np
from scipy import stats

DATA = _os.path.join(_REPO, "data")
RES = _os.path.join(_REPO, "results")
T_CRIT, N = 2.776, 5
HIGHER = {"TCR_pMHC_ipTM", "interface_pLDDT"}
LOWER = {"dG_kcal_mol", "Kd"}
ANIMAL = {"RHGEVNDSTTVEPVL", "GPGEVNDSTTVEPIL", "GPGEVNDSTTVEPVL"}

rows = list(csv.DictReader(open(_os.path.join(DATA, "per_seed_metrics.csv"))))
g = defaultdict(lambda: defaultdict(dict)); role = {}
for r in rows:
    if r["metric"] in HIGHER | LOWER and r["provenance"] == "WORKBOOK":
        k = (r["system"], r["peptide_label"], r["sequence"])
        g[k][r["metric"]][int(r["seed"])] = float(r["value"]); role[k] = r["role"]

def sv(k, m): return [g[k][m][s] for s in sorted(g[k][m])]
def mean(v): return sum(v)/len(v)
def sd(v):
    m = mean(v); return math.sqrt(sum((x-m)**2 for x in v)/(len(v)-1))

native_ci = {}
for k in g:
    if role[k] != "native_HIP": continue
    for m in HIGHER | LOWER:
        v = sv(k, m); mu = mean(v); hw = T_CRIT*sd(v)/math.sqrt(N)
        native_ci[(k[0], m)] = (mu-hw, mu+hw)

cands = [k for k in g if role[k] == "candidate"]

def gate(k, metrics):
    for m in metrics:
        mu = mean(sv(k, m)); lo, hi = native_ci[(k[0], m)]
        if m in HIGHER:
            if mu < lo: return False
        else:
            if mu > hi: return False
    return True

SETS = {
    "4-metric (as published: ipTM, i-pLDDT, dG, Kd)": ["TCR_pMHC_ipTM","interface_pLDDT","dG_kcal_mol","Kd"],
    "3-metric (drop Kd; dG counted once)":            ["TCR_pMHC_ipTM","interface_pLDDT","dG_kcal_mol"],
    "2-metric (ipTM + dG; one confidence, one energy)":["TCR_pMHC_ipTM","dG_kcal_mol"],
    "2-metric (i-pLDDT + dG)":                        ["interface_pLDDT","dG_kcal_mol"],
}
print("="*76); print("PART 1 -- does the gate's dimensionality change the nominated set?"); print("="*76)
results = {}
for name, ms in SETS.items():
    p = {k for k in cands if gate(k, ms)}
    nom = {k for k in p if k[2] not in ANIMAL}
    results[name] = nom
    print(f"{name:52s} pass={len(p):2d}  nominees={len(nom):2d}")
base = results["4-metric (as published: ipTM, i-pLDDT, dG, Kd)"]
print()
for name, nom in results.items():
    same = nom == base
    print(f"  identical to published set? {str(same):5s}  <- {name}")
    if not same:
        for d in sorted(base - nom): print(f"       lost: {d[2]} ({d[0]})")
        for d in sorted(nom - base): print(f"      added: {d[2]} ({d[0]})")

# ---------------------------------------------------------------- PART 2
print()
print("="*76); print("PART 2 -- audit of stored workbook values vs recomputation"); print("="*76)
summ = list(csv.DictReader(open(_os.path.join(DATA, "per_peptide_summary.csv"))))
ctrl = {k[0]: sv(k, "TCR_pMHC_ipTM") for k in g if role[k] == "negative_control"}

disc = []
for r in summ:
    k = (r["system"], r["peptide_label"], r["sequence"])
    if k not in g: continue
    iptm = sv(k, "TCR_pMHC_ipTM")

    # stored SD column
    try:
        stored_sd = float(r["tcr_pmhc_iptm_sd_stored"])
        real_sd = sd(iptm)
        if abs(stored_sd - real_sd) > 0.005:
            disc.append((k, "SD", stored_sd, round(real_sd,4),
                         "equals the MEAN" if abs(stored_sd-mean(iptm))<0.005 else ""))
    except (ValueError, KeyError): pass

    # stored p-value (candidates only)
    if role[k] == "candidate":
        try:
            sp = float(r["p_value_stored"])
            rp = stats.ttest_ind(iptm, ctrl[k[0]], equal_var=False).pvalue
            if abs(sp - rp) > 0.005:
                disc.append((k, "p-value", sp, round(float(rp),4), ""))
        except (ValueError, KeyError): pass

    # stored CI text "0.85 ± 0.01"
    txt = r.get("tcr_pmhc_iptm_ci_text","")
    if "±" in txt:
        try:
            m_s, h_s = [float(x.strip()) for x in txt.split("±")]
            hw = T_CRIT*sd(iptm)/math.sqrt(N)
            if abs(m_s-mean(iptm))>0.006 or abs(h_s-hw)>0.006:
                disc.append((k,"CI text",txt,f"{mean(iptm):.2f} ± {hw:.2f}",""))
        except ValueError: pass

by = defaultdict(list)
for d in disc: by[d[1]].append(d)
print(f"peptides audited: {len(summ)}   discrepancies: {len(disc)}\n")
for kind, items in sorted(by.items()):
    print(f"--- {kind}: {len(items)} ---")
    for k, _, stored, recomp, note in items[:25]:
        print(f"   {k[0]:22s} {k[2]:16s} stored={stored}  recomputed={recomp}  {note}")
    print()

with open(_os.path.join(RES,"workbook_audit.csv"),"w",newline="") as f:
    w=csv.writer(f); w.writerow(["system","peptide","sequence","field","stored","recomputed","note"])
    for k,kind,s_,r_,n_ in disc: w.writerow([k[0],k[1].strip(),k[2],kind,s_,r_,n_])
print(f"Wrote {RES}/workbook_audit.csv")
