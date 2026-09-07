#!/usr/bin/env python3
"""
Track B revision -- shuffled-decoy negatives (B1-B4).

B1: keep the 5 documented negatives as a separate gold-standard subset.
B2: generate shuffled decoys from the 36 cognate pairs -- each TCR paired with
    peptides from a DIFFERENT, unrelated source antigen, MHC allele held
    constant, excluding same-family and sequence-similar pairings. Target 100+.
B3: seed-row mislabeling (Cerosaletti 2024) already corrected in the delivered
    CSV (all six Clone_81/Clone_566/P196-1 rows are cognate) -- carried forward
    unchanged, noted here for the provenance record.
B4: allele composition reported (already in tcr_pmhc_benchmark_notes.md
    Section 3); this script does not repeat that analysis.

Also fixes a real data bug found while building this: row pair_id=11 (W321,
cognate) had the curation agent's caveat text concatenated into the
peptide_seq field itself, not just confidence_note. Corrected here, not
silently -- logged below.
"""
import csv
from itertools import combinations
import os as _os
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))

DATA = _os.path.join(_REPO, "data")

rows = list(csv.DictReader(open(f"{DATA}/tcr_pmhc_benchmark.csv")))

# ---- fix the W321 peptide_seq data bug ----
fixed = 0
for r in rows:
    if r["pair_id"] == "11" and r["tcr_source_id"] == "W321" and r["class"] == "cognate":
        assert r["peptide_seq"].startswith("QPFPQPEQPFP"), "unexpected content, aborting fix"
        old = r["peptide_seq"]
        r["peptide_seq"] = "QPFPQPEQPFP"
        r["confidence_note"] = (r["confidence_note"] + " | DATA FIX: peptide_seq field originally "
                                 "contained the curation caveat text concatenated after the real "
                                 "11-residue sequence; corrected here to the sequence alone. "
                                 f"Original raw value: {old!r}")
        fixed += 1
print(f"Fixed {fixed} row(s) with corrupted peptide_seq (expected 1: W321).")

cognate = [r for r in rows if r["class"] == "cognate"]
documented_negatives = [r for r in rows if r["class"] == "non_cognate"]
print(f"Cognate rows: {len(cognate)}  |  Documented negatives (B1, kept separate): {len(documented_negatives)}")

# ---- hand-curated antigen-family tags (documented, not automated/guessed) ----
# Keyed by (tcr_source_id, peptide_name) since a few names repeat across rows.
FAMILY = {
    "GAD65_377-396": "GAD65",
    "MP54_97-116": "Influenza_M1",
    "DQ2.5-glia-alpha1a": "Gliadin",
    "DQ2.5-glia-omega1": "Gliadin",
    "DQ2.5-glia-alpha2 (deamidated)": "Gliadin",
    "DQ2.5-glia-alpha2 (non-deamidated)": "Gliadin",
    "P.fluorescens SGDS-derived bacterial mimic of DQ2.5-glia-alpha1a": "Gliadin",  # explicit molecular mimic
    "Influenza HA1 306-318 (HA1.7 epitope)": "Influenza_HA",
    "Influenza HA1 306-318 (HA 1.7 epitope)": "Influenza_HA",
    "DQ8-glia-alpha1 (gliadin, MM1 peptide)": "Gliadin",
    "Insulin B:9-23 mimotope p8E9E (8E9E11ss engineered superagonist)": "Insulin",
    "Chicken conalbumin (CA) peptide": "Conalbumin",
    "Myelin basic protein Ac1-11 (EAE autoantigen)": "MBP",
    "Hybrid insulin peptide (HIP, chromogranin A/insulin C-peptide fusion)": "Insulin_HIP",
    "Insulin B-chain mimotope, register-3B (Reg3B, p8G9E)": "Insulin",
    "Mutant triosephosphate isomerase (mutTPI) tumor neoantigen": "TPI_neoantigen",
    "DQ2.5-glia-alpha1a (deamidated)": "Gliadin",
    "Proinsulin C-peptide PI40-54 (register G9EL11C stabilized)": "Insulin_HIP",
    "Hybrid insulin peptide HIP3 (PI40-47/NPY68-74, L11C stabilized)": "Insulin_HIP",
    "Hybrid insulin peptide HIP1 (PI40-47/IAPP74-80, L11C stabilized)": "Insulin_HIP",
    "Citrullinated vimentin 64cit(59-71)": "Cit_Vimentin",
    "Citrullinated alpha-enolase 15cit(10-22)": "Cit_Enolase",
    "Citrullinated tenascin-C TNC1014,1016cit (1013-1024, double citrulline)": "Cit_TenascinC",
    "Myelin basic protein MBP85-99": "MBP",
    "Myelin oligodendrocyte glycoprotein MOG35-55 (mouse EAE)": "MOG",
    "Chromogranin-A fragment (WE14-related, mouse T1D)": "Insulin_HIP",
    "Thyroid peroxidase self-peptide (mouse Graves' disease model)": "TPO",
}

missing_family = [r["peptide_name"] for r in cognate if r["peptide_name"] not in FAMILY]
if missing_family:
    print("WARNING -- peptides with no family tag (would be excluded from decoy generation):")
    for m in set(missing_family):
        print(" ", m)

def strip_mods(seq):
    return seq.replace("(cit)", "X")

def shares_5mer(a, b):
    a, b = strip_mods(a), strip_mods(b)
    if len(a) < 5 or len(b) < 5:
        return len(a) >= 3 and len(b) >= 3 and a in b or b in a
    a_kmers = {a[i:i+5] for i in range(len(a) - 4)}
    b_kmers = {b[i:i+5] for i in range(len(b) - 4)}
    return len(a_kmers & b_kmers) > 0

# ---- build decoys ----
decoys = []
for tcr_row in cognate:
    tcr_id = tcr_row["tcr_source_id"]
    tcr_cdr3a, tcr_cdr3b = tcr_row["CDR3a"], tcr_row["CDR3b"]
    mhc = tcr_row["MHC_beta_allele"]
    own_family = FAMILY.get(tcr_row["peptide_name"])
    own_peptides = {strip_mods(r["peptide_seq"]) for r in cognate if r["tcr_source_id"] == tcr_id}
    own_families = {FAMILY.get(r["peptide_name"]) for r in cognate if r["tcr_source_id"] == tcr_id}

    for pep_row in cognate:
        if pep_row["tcr_source_id"] == tcr_id:
            continue  # never pair a TCR with its own peptide set
        pep_family = FAMILY.get(pep_row["peptide_name"])
        if pep_family is None or pep_family in own_families:
            continue  # unclassified or same/related antigen family -> excluded
        if pep_row["MHC_beta_allele"] != mhc:
            continue  # MHC allele must match
        pep_seq = strip_mods(pep_row["peptide_seq"])
        if any(shares_5mer(pep_seq, op) for op in own_peptides):
            continue  # excluded: shares a 5-mer with one of the TCR's own cognate peptides
        decoys.append(dict(
            pair_id=f"decoy_{len(decoys)+1}",
            class_="shuffled_decoy",
            tcr_source_id=tcr_id,
            TRAV=tcr_row["TRAV"], TRAJ=tcr_row["TRAJ"], CDR3a=tcr_cdr3a,
            TRBV=tcr_row["TRBV"], TRBJ=tcr_row["TRBJ"], CDR3b=tcr_cdr3b,
            peptide_name=pep_row["peptide_name"], peptide_seq=pep_row["peptide_seq"],
            MHC_class="II", MHC_alpha_allele=tcr_row["MHC_alpha_allele"],
            MHC_beta_allele=mhc,
            MHC_alpha_seq=tcr_row["MHC_alpha_seq"], MHC_beta_seq=tcr_row["MHC_beta_seq"],
            evidence_type="shuffled_decoy_assumed_negative",
            source_id=f"combinatorial: TCR from {tcr_row['source_id']}; peptide from {pep_row['source_id']}",
            source_url="",
            confidence_note=(
                f"ASSUMED negative, not experimentally verified. TCR {tcr_id} (cognate for "
                f"{own_family} family, MHC {mhc}) paired with unrelated-family peptide "
                f"'{pep_row['peptide_name']}' ({pep_family} family), same MHC allele, no shared "
                "5-mer with any of this TCR's own cognate peptides. Based on the assumption that "
                "TCR specificity is narrow enough that an unrelated-antigen, sequence-dissimilar "
                "peptide is very unlikely to be a true binder -- this is a documented assumption, "
                "not a verified fact; false-negative rate in this set is expected to be small but "
                "is not zero and has not been measured."
            ),
        ))

print(f"\nShuffled decoys generated: {len(decoys)}")

fieldnames = list(cognate[0].keys())
with open(f"{DATA}/tcr_pmhc_benchmark_decoys.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for d in decoys:
        row = {k: d.get(k, d.get(k.replace("class", "class_"), "")) for k in fieldnames}
        row["class"] = "shuffled_decoy"
        w.writerow(row)

# also write the corrected main benchmark CSV (with W321 fix applied) back out
with open(f"{DATA}/tcr_pmhc_benchmark.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print(f"Wrote {DATA}/tcr_pmhc_benchmark_decoys.csv ({len(decoys)} rows)")
print(f"Wrote corrected {DATA}/tcr_pmhc_benchmark.csv (W321 peptide_seq fixed)")

# summary by TCR
from collections import Counter
per_tcr = Counter(d["tcr_source_id"] for d in decoys)
print("\nDecoys per source TCR:")
for tcr, n in sorted(per_tcr.items()):
    print(f"  {tcr:12s} {n}")
