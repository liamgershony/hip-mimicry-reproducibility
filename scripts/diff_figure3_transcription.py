#!/usr/bin/env python3
"""
Phase 1.3 -- full cell-by-cell diff of the OCR/screenshot-transcribed data used
to build the delivered Figure 3 (per-system ipTM distribution) against the real
Master Results Sheet.xlsx, via per_seed_metrics.csv.

Every discrepancy is reported. Nothing is silently corrected here.
"""
import csv
from collections import defaultdict
import os as _os
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))

OLD_SCRIPT = "/private/tmp/claude-501/-Users-liamgershony-Downloads-Does-it-Mimic/827d0b88-9cbc-4225-ac6c-c7b11a662013/scratchpad/build_v3_figures.py"
DATA = _os.path.join(_REPO, "data")
RESULTS = _os.path.join(_REPO, "results")

# --- extract just the data-definition portion of the old script (lines before
# the SYSTEMS assignment, i.e. before any plotting code runs) ---
with open(OLD_SCRIPT) as f:
    lines = f.readlines()
data_src = "".join(lines[:106])  # through the SYSTEMS = [...] line
ns = {}
exec(data_src, ns)
HIP6, HIP11_cis, HIP11_trans = ns["HIP6"], ns["HIP11_cis"], ns["HIP11_trans"]

OLD_BLOCKS = [
    ("HIP6/A2.11/DQ8", HIP6),
    ("HIP11/8.E3/DQ8", HIP11_cis),
    ("HIP11/8.E3/DQ8-trans", HIP11_trans),
]

# --- load the real data, in the same positional candidate order (row order in
# the workbook, which is how both the old transcription and the new CSV were
# both ultimately ordered -- confirmed identical sequence order in Phase 0) ---
rows = list(csv.DictReader(open(f"{DATA}/per_seed_metrics.csv")))
iptm_rows = [r for r in rows if r["metric"] == "TCR_pMHC_ipTM"]
grouped = defaultdict(dict)
order = defaultdict(list)
for r in iptm_rows:
    key = (r["system"], r["peptide_label"], r["sequence"])
    grouped[key][int(r["seed"])] = float(r["value"])

# real candidates in workbook row order, per system
real_by_system = defaultdict(list)
seen = defaultdict(set)
for r in iptm_rows:
    if r["role"] != "candidate":
        continue
    key = (r["system"], r["peptide_label"], r["sequence"])
    if key not in seen[r["system"]]:
        seen[r["system"]].add(key)
        real_by_system[r["system"]].append(key)

discrepancies = []
total_compared = 0
for system, old_block in OLD_BLOCKS:
    old_cands = old_block["candidates"]  # list of (iptm5, rmsd, pass_bool), in row order
    real_cands = real_by_system[system]
    if len(old_cands) != len(real_cands):
        print(f"WARNING: candidate count mismatch for {system}: old={len(old_cands)} real={len(real_cands)}")
    for i, ((old_iptm5, old_rmsd, old_pass), real_key) in enumerate(zip(old_cands, real_cands)):
        real_vals = [grouped[real_key][s] for s in range(1, 6)]
        for seed_idx in range(5):
            total_compared += 1
            old_v = old_iptm5[seed_idx]
            real_v = real_vals[seed_idx]
            if abs(old_v - real_v) > 1e-9:
                discrepancies.append(dict(
                    system=system, candidate_index=i + 1, sequence=real_key[2],
                    seed=seed_idx + 1, transcribed_value=old_v, workbook_value=real_v,
                    diff=round(real_v - old_v, 4),
                ))

with open(f"{RESULTS}/figure3_transcription_diff.csv", "w", newline="") as f:
    if discrepancies:
        w = csv.DictWriter(f, fieldnames=list(discrepancies[0].keys()))
        w.writeheader()
        w.writerows(discrepancies)

print(f"Total ipTM values compared: {total_compared}")
print(f"Discrepancies found: {len(discrepancies)}")
print()
for d in discrepancies:
    print(f"  {d['system']:22s} cand#{d['candidate_index']:2d} {d['sequence']:16s} seed{d['seed']}: "
          f"transcribed={d['transcribed_value']}  workbook={d['workbook_value']}  diff={d['diff']:+.2f}")
