#!/usr/bin/env python3
"""
Build the TCRmodel2 submission sheet: one row per JOB to run, with empty
columns for every metric to be filled in, plus a job-URL column.

Inputs (all already on disk, nothing invented):
  data/per_seed_metrics.csv          -- the 61 manuscript systems + published values
  data/tcr_pmhc_benchmark.csv        -- 36 cognate + 5 documented negatives
  data/tcr_pmhc_benchmark_decoys.csv -- 100 shuffled decoys

TCR clone -> gene/CDR3 assignments for the manuscript's own systems come from
Table 1 of the manuscript, and were independently verified in Phase 0 against
the real AlphaFold 3 job CIF files (entity_poly sequences matched byte-for-byte
against TCRmodel2's own TRAV/TRAJ/TRBV/TRBJ reference sequences). MHC allele
assignments likewise verified against the CIF entity 3/4 sequences.
"""
import csv
import math
from collections import defaultdict
import os as _os
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))

DATA = _os.path.join(_REPO, "data")
OUT = f"{DATA}/tcrmodel2_submission_sheet.csv"

# ---------------------------------------------------------------------------
# Verified TCR + MHC definitions for the manuscript's three candidate-bearing
# systems. Sequences are TCRmodel2's own dropdown values (allele -> sequence),
# confirmed identical to the deposited AF3 job CIF entity sequences.
# ---------------------------------------------------------------------------
DQA_0301 = "IVADHVASYGVNLYQSYGPSGQYSHEFDGDEEFYVDLERKETVWQLPLFRRFRRFDPQFALTNIAVLKHNLNIVIKRSNSTAAT"
DQA_0501 = "IVADHVASYGVNLYQSYGPSGQYTHEFDGDEQFYVDLGRKETVWCLPVLRQFRFDPQFALTNIAVLKHNLNSLIKRSNSTAAT"
DQB_0302 = "SPEDFVYQFKGMCYFTNGTERVRLVTRYIYNREEYARFDSDVGVYRAVTPLGPPAAEYWNSQKEVLERTRAELDTVCRHNYQLELRTTLQ"

SYSTEMS = {
    "HIP6/A2.11/DQ8": dict(
        clone="A2.11",
        TRAV="TRAV38-1*01", TRAJ="TRAJ54*01", CDR3a="CAFMGAGAQKLVF",
        TRBV="TRBV4-3*01", TRBJ="TRBJ2-3*01", CDR3b="CASSQILRGGPPDTQYF",
        mhc_a_allele="HLA-DQA1*03:01", mhc_b_allele="HLA-DQB1*03:02",
        mhc_a_seq=DQA_0301, mhc_b_seq=DQB_0302,
        restriction="DQ8",
    ),
    "HIP11/8.E3/DQ8": dict(
        clone="8.E3",
        TRAV="TRAV2*01", TRAJ="TRAJ37*01", CDR3a="CAVDGSGNTGKLIF",
        TRBV="TRBV4-1*01", TRBJ="TRBJ2-7*01", CDR3b="CASSQDLAGVREQYF",
        mhc_a_allele="HLA-DQA1*03:01", mhc_b_allele="HLA-DQB1*03:02",
        mhc_a_seq=DQA_0301, mhc_b_seq=DQB_0302,
        restriction="DQ8 (cis)",
    ),
    "HIP11/8.E3/DQ8-trans": dict(
        clone="8.E3",
        TRAV="TRAV2*01", TRAJ="TRAJ37*01", CDR3a="CAVDGSGNTGKLIF",
        TRBV="TRBV4-1*01", TRBJ="TRBJ2-7*01", CDR3b="CASSQDLAGVREQYF",
        mhc_a_allele="HLA-DQA1*05:01", mhc_b_allele="HLA-DQB1*03:02",
        mhc_a_seq=DQA_0501, mhc_b_seq=DQB_0302,
        restriction="DQ8-trans",
    ),
}

# Pre-declared pilot set (Track C2), by (system, sequence)
PILOT = {
    ("HIP11/8.E3/DQ8-trans", "RRNVATLQAENVTG"): "high-scoring nominated candidate (lead, SARS-CoV-2)",
    ("HIP11/8.E3/DQ8-trans", "REDTVSVKSEPVSE"): "mid-scoring candidate (passed Welch, failed CI gate)",
    ("HIP6/A2.11/DQ8", "GVEDVYLAGALEAQ"): "gate-failed candidate (failed Welch, p=0.4072)",
    ("HIP11/8.E3/DQ8-trans", "SLQPLALEAEDLQV"): "native HIP11 reference",
    ("HIP11/8.E3/DQ8-trans", "KDVDAAVDAEVVQF"): "negative control (EBV AG876)",
}

# Tolerance bands from scripts/derive_tolerance_band.py (median / conservative)
BAND_IPTM_MED, BAND_IPTM_CONS = 0.036, 0.369
BAND_IPLDDT_MED, BAND_IPLDDT_CONS = 2.65, 15.71

# ---------------------------------------------------------------------------
# Load published per-seed values for the 61 manuscript systems
# ---------------------------------------------------------------------------
rows = list(csv.DictReader(open(f"{DATA}/per_seed_metrics.csv")))
pub = defaultdict(dict)
meta = {}
for r in rows:
    if r["provenance"] != "WORKBOOK":
        continue
    key = (r["system"], r["peptide_label"], r["sequence"])
    pub[key].setdefault(r["metric"], {})[int(r["seed"])] = float(r["value"])
    meta[key] = dict(role=r["role"], organism=r["source_organism"])

def mean_sd(vals):
    n = len(vals)
    m = sum(vals) / n
    sd = math.sqrt(sum((v - m) ** 2 for v in vals) / (n - 1)) if n > 1 else 0.0
    return m, sd

# ---------------------------------------------------------------------------
# Column layout
# ---------------------------------------------------------------------------
INPUT_COLS = [
    "job_id", "batch", "run_order", "submission_ready", "blocker",
    "system", "role", "peptide_label", "source_organism", "pilot_rationale",
    "tcr_clone", "TRAV", "TRAJ", "CDR3a", "TRBV", "TRBJ", "CDR3b",
    "peptide_seq", "mhc_class", "mhc_species",
    "mhc_alpha_allele", "mhc_beta_allele", "mhc_alpha_seq", "mhc_beta_seq",
    "amber_relax", "n_models_expected",
]
PUBLISHED_COLS = [
    "published_mean_tcr_pmhc_iptm", "published_sd_tcr_pmhc_iptm",
    "published_mean_interface_plddt", "published_sd_interface_plddt",
    "published_mean_dG_kcal_mol", "published_mean_Kd_M",
    "tolerance_band_iptm", "tolerance_band_interface_plddt", "tolerance_band_basis",
]
TRACKING_COLS = [
    "job_url", "job_id_returned_by_server", "server_url", "date_submitted",
    "tcrmodel2_version_string", "alphafold_backend_version", "random_seed_or_range",
    "wallclock_minutes", "run_notes",
]
def metric_cols():
    # Column names match TCRmodel2's own statistics.json keys exactly, so the
    # server output can be mapped straight across. Models are ranked_0..ranked_4
    # as emitted by the pipeline (NOT 1-indexed).
    cols = []
    for m in range(5):
        cols += [
            f"ranked{m}_tcr_pmhc_iptm",      # 'tcr-pmhc_iptm' -- the paper's primary metric
            f"ranked{m}_interface_plddt",    # computed by calc_iplddt, chn1=ABC vs chn2=DE
            f"ranked{m}_iptm",               # overall iptm
            f"ranked{m}_ptm",
            f"ranked{m}_plddt",
            f"ranked{m}_ranking_confidence", # 'model confidence'
        ]
    return cols
RMSD_COLS = ["rmsd_angstrom", "rmsd_model_used"]

FIELDNAMES = INPUT_COLS + PUBLISHED_COLS + TRACKING_COLS + metric_cols() + RMSD_COLS

def blank_row():
    return {k: "" for k in FIELDNAMES}

out_rows = []
job_n = 0

# ---------------------------------------------------------------------------
# Batch 1 + 2: the 61 manuscript systems
# ---------------------------------------------------------------------------
manuscript_jobs = []
for key, metrics in pub.items():
    system, label, seq = key
    if system not in SYSTEMS:
        continue
    iptm = [metrics["TCR_pMHC_ipTM"][s] for s in sorted(metrics["TCR_pMHC_ipTM"])]
    ipl = [metrics["interface_pLDDT"][s] for s in sorted(metrics["interface_pLDDT"])]
    dg = [metrics["dG_kcal_mol"][s] for s in sorted(metrics["dG_kcal_mol"])]
    kd = [metrics["Kd"][s] for s in sorted(metrics["Kd"])]
    manuscript_jobs.append((key, iptm, ipl, dg, kd))

# pilot first, then the rest in system/label order
def sort_key(item):
    (system, label, seq), *_ = item
    is_pilot = (system, seq) in PILOT
    try:
        n = int("".join(c for c in label if c.isdigit()) or 0)
    except ValueError:
        n = 0
    return (0 if is_pilot else 1, system, n, label)

for (key, iptm, ipl, dg, kd) in sorted(manuscript_jobs, key=sort_key):
    system, label, seq = key
    S = SYSTEMS[system]
    is_pilot = (system, seq) in PILOT
    job_n += 1
    r = blank_row()
    r.update(
        job_id=f"J{job_n:03d}",
        batch="1_pilot" if is_pilot else "2_tier1_regeneration",
        run_order=job_n,
        submission_ready="yes",
        blocker="",
        system=system,
        role=meta[key]["role"],
        peptide_label=label.strip(),
        source_organism=meta[key]["organism"],
        pilot_rationale=PILOT.get((system, seq), ""),
        tcr_clone=S["clone"],
        TRAV=S["TRAV"], TRAJ=S["TRAJ"], CDR3a=S["CDR3a"],
        TRBV=S["TRBV"], TRBJ=S["TRBJ"], CDR3b=S["CDR3b"],
        peptide_seq=seq,
        mhc_class="II", mhc_species="human",
        mhc_alpha_allele=S["mhc_a_allele"], mhc_beta_allele=S["mhc_b_allele"],
        mhc_alpha_seq=S["mhc_a_seq"], mhc_beta_seq=S["mhc_b_seq"],
        amber_relax="TRUE", n_models_expected=5,
    )
    m_iptm, sd_iptm = mean_sd(iptm)
    m_ipl, sd_ipl = mean_sd(ipl)
    m_dg, _ = mean_sd(dg)
    m_kd, _ = mean_sd(kd)
    high_var = sd_iptm > 0.10
    r.update(
        published_mean_tcr_pmhc_iptm=round(m_iptm, 4),
        published_sd_tcr_pmhc_iptm=round(sd_iptm, 4),
        published_mean_interface_plddt=round(m_ipl, 3),
        published_sd_interface_plddt=round(sd_ipl, 3),
        published_mean_dG_kcal_mol=round(m_dg, 3),
        published_mean_Kd_M=f"{m_kd:.3e}",
        tolerance_band_iptm=BAND_IPTM_CONS if high_var else BAND_IPTM_MED,
        tolerance_band_interface_plddt=BAND_IPLDDT_CONS if high_var else BAND_IPLDDT_MED,
        tolerance_band_basis="conservative (published per-seed SD > 0.10)" if high_var
                              else "median-case",
        server_url="https://tcrmodel.ibbr.umd.edu/",
    )
    out_rows.append(r)

# ---------------------------------------------------------------------------
# Batches 3-5: the Phase 2 benchmark
# ---------------------------------------------------------------------------
def add_benchmark(path, batch_name):
    global job_n
    for b in csv.DictReader(open(path)):
        job_n += 1
        r = blank_row()
        genes = [b["TRAV"], b["TRAJ"], b["TRBV"], b["TRBJ"]]
        missing_genes = [g for g in genes if not g or g.upper().startswith("NOT FOUND")]
        no_mhc_seq = (not b["MHC_alpha_seq"] or b["MHC_alpha_seq"].upper().startswith("NOT FOUND")
                      or not b["MHC_beta_seq"] or b["MHC_beta_seq"].upper().startswith("NOT FOUND"))
        is_mouse = "mouse" in b["MHC_beta_allele"].lower() or b["MHC_beta_allele"].startswith("I-A")
        blockers = []
        if missing_genes:
            blockers.append(f"{len(missing_genes)} of 4 V/J gene names NOT FOUND (CDR3 verified) "
                            "-- submit via full-chain 'Enter sequence' route or resolve genes first")
        if no_mhc_seq:
            blockers.append("MHC chain sequence(s) NOT FOUND -- need IMGT/HLA or UniProt lookup")
        if is_mouse:
            blockers.append("mouse MHC II -- confirm allele is available in TCRmodel2's mouse dropdown")
        r.update(
            job_id=f"J{job_n:03d}",
            batch=batch_name,
            run_order=job_n,
            submission_ready="no" if blockers else "yes",
            blocker=" | ".join(blockers),
            system=f"benchmark:{b['tcr_source_id']}",
            role=b["class"],
            peptide_label=b["peptide_name"],
            source_organism=b.get("evidence_type", ""),
            tcr_clone=b["tcr_source_id"],
            TRAV=b["TRAV"], TRAJ=b["TRAJ"], CDR3a=b["CDR3a"],
            TRBV=b["TRBV"], TRBJ=b["TRBJ"], CDR3b=b["CDR3b"],
            peptide_seq=b["peptide_seq"],
            mhc_class="II",
            mhc_species="mouse" if is_mouse else "human",
            mhc_alpha_allele=b["MHC_alpha_allele"], mhc_beta_allele=b["MHC_beta_allele"],
            mhc_alpha_seq=b["MHC_alpha_seq"], mhc_beta_seq=b["MHC_beta_seq"],
            amber_relax="TRUE", n_models_expected=5,
            server_url="https://tcrmodel.ibbr.umd.edu/",
            run_notes=b.get("source_id", "")[:200],
        )
        out_rows.append(r)

# split the main benchmark file by class so batches are clean
cog_path = f"{DATA}/_tmp_cognate.csv"
neg_path = f"{DATA}/_tmp_negative.csv"
bench = list(csv.DictReader(open(f"{DATA}/tcr_pmhc_benchmark.csv")))
bf = list(bench[0].keys())
for path, cls in [(cog_path, "cognate"), (neg_path, "non_cognate")]:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=bf)
        w.writeheader()
        w.writerows([b for b in bench if b["class"] == cls])

add_benchmark(cog_path, "3_benchmark_cognate")
add_benchmark(neg_path, "4_benchmark_documented_negative")
add_benchmark(f"{DATA}/tcr_pmhc_benchmark_decoys.csv", "5_benchmark_shuffled_decoy")

import os
os.remove(cog_path); os.remove(neg_path)

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=FIELDNAMES)
    w.writeheader()
    w.writerows(out_rows)

# ---------------------------------------------------------------------------
from collections import Counter
print(f"Wrote {OUT}")
print(f"Total jobs: {len(out_rows)}   Total columns: {len(FIELDNAMES)}")
print()
for b, n in sorted(Counter(r["batch"] for r in out_rows).items()):
    ready = sum(1 for r in out_rows if r["batch"] == b and r["submission_ready"] == "yes")
    print(f"  {b:35s} {n:4d} jobs   ({ready} submission-ready, {n-ready} blocked)")
print()
print(f"Submission-ready overall: {sum(1 for r in out_rows if r['submission_ready']=='yes')} / {len(out_rows)}")
print()
print("Pilot jobs:")
for r in out_rows:
    if r["batch"] == "1_pilot":
        print(f"  {r['job_id']}  {r['peptide_seq']:16s} {r['system']:22s} "
              f"pub ipTM={r['published_mean_tcr_pmhc_iptm']}  band=±{r['tolerance_band_iptm']}")
