# Track B revision — shuffled-decoy negatives (B1–B4)

## B1. Documented negatives kept as gold-standard subset

The original 5 documented negatives remain in `tcr_pmhc_benchmark.csv`, `class=non_cognate`,
unchanged. Report these separately from the decoys below — they are the only rows in this
whole benchmark with an actual experimental negative result behind them.

## B2. Shuffled decoys — `data/tcr_pmhc_benchmark_decoys.csv`, 100 rows, `class=shuffled_decoy`

**Generation rule** (`scripts/build_shuffled_decoys.py`): for each of the 36 cognate TCRs,
paired with every peptide from a *different* cognate row that simultaneously satisfies all
three:
1. **Same MHC β-allele** as the TCR's own restriction (a decoy across MHC alleles would
   confound TCR specificity with MHC-binding, not test what we want).
2. **Different, unrelated source-antigen family** — hand-curated family tags for all 27
   distinct cognate peptides (table in the script; e.g. all gliadin-derived epitopes,
   including one paper's explicit bacterial molecular-mimic of a gliadin epitope, are one
   "Gliadin" family; all HIP/proinsulin/chromogranin-related peptides are one "Insulin_HIP"
   family — this one matters most, since it's adjacent to this manuscript's own subject).
3. **No shared 5-residue exact substring** with any of that TCR's own cognate peptides
   (citrullination markers normalized before comparison) — a second, sequence-level check
   independent of the family tags, catching near-duplicates the family curation might miss.

**Result: exactly 100 decoys**, hitting the target. Coverage is uneven by design, not by
mistake — human DRB1\*04:01 and DQ2.5-restricted TCRs (Clone_81, Clone_566, P196-1, HA1.7,
A03, A07, RA2.7, PB) get 9–11 decoys each, because those alleles have many candidate peptides
in the cognate set to draw from. The four mouse-MHC TCRs (4.1, 8F10, BDC-2.5, D10, G4, TAZ10)
get only 1–2 each, because their alleles (I-Ag7, I-Au, I-Ak, I-Ab) have very few other entries
to pair against. **This decoy set is not allele-balanced** — flagging now, consistent with B4.

**The stated assumption, not a verified fact**: these 100 pairings are *assumed* negative on
the basis that TCR specificity is narrow enough that an unrelated-antigen, sequence-dissimilar
peptide on the correct MHC is very unlikely to be a true binder. This is a reasonable
assumption grounded in what's known about TCR degeneracy, but it is not itself measured here
— the false-negative rate of this decoy set (i.e., how many of the 100 might actually be
cryptic true binders) is unknown. Every decoy row's `confidence_note` states this explicitly
rather than implying decoy status carries the same evidence weight as the 5 documented
negatives or the 36 cognate positives.

## B3. Seed-row correction (carried forward, not redone)

Confirmed already correctly applied in the delivered CSV: all six Clone_81/Clone_566/P196-1
rows are `cognate`, per Cerosaletti et al. 2024 (PMID 38871688) showing cross-reactivity
rather than one cognate/one non-cognate peptide per clone. No seed row is coded
`non_cognate` or `shuffled_decoy`.

## B4. Allele composition — for Phase 2 to act on, not resolved here

From the existing `tcr_pmhc_benchmark_notes.md` §3: DRB1\*04:01 and DQ2.5 together are 59%
of the 41 cognate+documented-negative rows. **This skew is arguably favorable for this
project specifically** — the manuscript's own systems are DQ8- and DQ2-restricted, and
DQ8/DQ2.5 share substantial peptide-binding-groove architecture with DQ2.5 being the single
best-represented class-II system in this benchmark. But favorable-for-us is not the same as
representative, and pooling across alleles would hide allele-specific effects either way.
**Requirement carried into Phase 2, not executed now**: report AUROC/PR stratified by MHC
allele (at minimum: DRB1\*04:01, DQ2.5, DQ8, and "all mouse alleles pooled" as a fourth
bucket given their small individual n) alongside the pooled number, not the pooled number
alone.

## Data bug found and fixed while building this (not part of B1–B4, but real)

`pair_id=11` (W321, cognate, DQ2.5-glia-omega1) had the curation agent's caveat text
concatenated directly into the `peptide_seq` field — the cell literally contained
`"QPFPQPEQPFP (registration per DQ2.5-glia-omega1 epitope; exact crystallized flanks not
available -- no structure solved for W321)"` instead of just the 11-residue sequence. Fixed
in `tcr_pmhc_benchmark.csv` (now `QPFPQPEQPFP`); the full original text and the fix are both
recorded in that row's `confidence_note` for the audit trail. This would have broken any
downstream modeling submission for that row had it gone unnoticed — worth double-checking
the rest of the CSV's peptide_seq fields by eye before any of these get submitted anywhere.
