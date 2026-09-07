#!/usr/bin/env python3
"""
Build the PRODIGY worklist: one row per MODEL PDB (not per TCRmodel2 job),
because PRODIGY is run separately on each structure.

202 TCRmodel2 jobs x 5 models = 1,010 PRODIGY runs.
For the 61 manuscript systems alone: 61 x 5 = 305, which is exactly the
"305 seed values" the manuscript's own Section 3.15 refers to.

CHAIN SELECTION -- the parameter PRODIGY actually needs, and the one thing that
would silently corrupt every dG if set wrong.

Derived by tracing tcr_utils.renumber_pdb() in piercelab/tcrmodel2, which runs
unconditionally on every output model (run_tcrmodel2.py "Renumber output"
section). For mhc_class == 2 the chain letters are permuted from the input
FASTA order (TCRa, TCRb, Peptide, MHCa, MHCb -> A,B,C,D,E) as follows:

    A(TCRa) -> F ;  B(TCRb) -> G ;  D(MHCa) -> A ;  E(MHCb) -> B
    F -> D       ;  G -> E

giving the FINAL output-PDB mapping:

    A = MHC alpha
    B = MHC beta
    C = peptide
    D = TCR alpha
    E = TCR beta

This is corroborated inside the same codebase: calc_iplddt() uses
chn1 = "ABC" (pMHC) vs chn2 = "DE" (TCR), and the CDR3 b-factor code selects
'chainID D' for CDR3a and 'chainID E' for CDR3b. It also matches the
manuscript's Section 2.7 wording ("the interface between the TCR alpha/beta
chains and the pMHC complex").

    => PRODIGY selection:  group 1 = D,E   (TCR)
                           group 2 = A,B,C (pMHC)

STILL WORTH SPOT-CHECKING against the first pilot PDB actually returned by the
web server, since the server build may differ from the public repo. Verifying
one file settles it for all 1,010.
"""
import csv
import os as _os
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))

DATA = _os.path.join(_REPO, "data")
SRC = f"{DATA}/tcrmodel2_submission_sheet.csv"
OUT = f"{DATA}/prodigy_worklist.csv"

# PRODIGY default temperature is 25 C. Phase 1's Kd/dG consistency check
# recovered a median implied RT of 0.59208 kcal/mol from the published data,
# which corresponds to ~298 K -- i.e. the original work used the default and
# we should too, or the Kd column will not be comparable.
TEMPERATURE_C = 25

FIELDNAMES = [
    "prodigy_run_id", "job_id", "batch", "model_rank", "model_pdb_filename",
    "system", "role", "peptide_label", "peptide_seq",
    "prodigy_selection_group1_tcr", "prodigy_selection_group2_pmhc",
    "temperature_C", "server_url",
    # published reference (manuscript systems only; per-seed order is arbitrary
    # so compare distributions, not element-wise -- see note in README)
    "published_mean_dG_kcal_mol", "published_mean_Kd_M",
    # ---- to be filled in ----
    "prodigy_dG_kcal_mol", "prodigy_Kd_M", "n_intermolecular_contacts",
    "prodigy_job_url", "prodigy_version", "date_run", "notes",
]

jobs = list(csv.DictReader(open(SRC)))
out = []
n = 0
for j in jobs:
    for rank in range(5):
        n += 1
        out.append({
            "prodigy_run_id": f"P{n:04d}",
            "job_id": j["job_id"],
            "batch": j["batch"],
            "model_rank": rank,
            "model_pdb_filename": f"ranked_{rank}.pdb",
            "system": j["system"],
            "role": j["role"],
            "peptide_label": j["peptide_label"],
            "peptide_seq": j["peptide_seq"],
            "prodigy_selection_group1_tcr": "D,E",
            "prodigy_selection_group2_pmhc": "A,B,C",
            "temperature_C": TEMPERATURE_C,
            "server_url": "https://wenmr.science.uu.nl/prodigy/",
            "published_mean_dG_kcal_mol": j["published_mean_dG_kcal_mol"],
            "published_mean_Kd_M": j["published_mean_Kd_M"],
            "prodigy_dG_kcal_mol": "", "prodigy_Kd_M": "",
            "n_intermolecular_contacts": "", "prodigy_job_url": "",
            "prodigy_version": "", "date_run": "", "notes": "",
        })

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=FIELDNAMES)
    w.writeheader()
    w.writerows(out)

from collections import Counter
print(f"Wrote {OUT}")
print(f"Total PRODIGY runs: {len(out)}  ({len(jobs)} jobs x 5 models)")
print()
for b, c in sorted(Counter(r["batch"] for r in out).items()):
    print(f"  {b:35s} {c:5d} PRODIGY runs")
print()
man = sum(c for b, c in Counter(r["batch"] for r in out).items()
          if b in ("1_pilot", "2_tier1_regeneration"))
print(f"Manuscript systems (batches 1+2): {man} runs "
      f"-- matches the paper's own '305 seed values' (Sec 3.15): {man == 305}")
print(f"Pilot alone: {sum(1 for r in out if r['batch']=='1_pilot')} runs")
