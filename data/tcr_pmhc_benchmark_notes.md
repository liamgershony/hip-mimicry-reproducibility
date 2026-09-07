# TCR–pMHC-II Benchmark: Curation Notes

## 1. Final counts (real numbers, not padded to target)

| | count |
|---|---|
| **Cognate pairs** | **36** |
| **Non-cognate pairs** | **5** |
| **Total rows** | **41** |
| Distinct TCR clones (`tcr_source_id`) | 29 |
| TCR clones with a linked cognate **and** non-cognate pair | 5 (`W316`, `W321`, `S2`, `A2.13`, `PB`) |

The cognate target (≥20) was exceeded. **The non-cognate target (≥20) was not met — only 5 non-cognate pairs could be built**, and this is reported honestly rather than padded. Every non-cognate row is a genuine, citable "same TCR + same MHC restriction + different peptide + explicit negative result" pairing; see §4 for why more could not be found.

Evidence-type breakdown: 24 rows are `solved_structure` (a deposited PDB ternary complex whose chain composition and CDR3 sequences I parsed myself directly from the RCSB-deposited FASTA), 17 are `functional_assay` (proliferation, cytokine release, SPR/tetramer binding, or CD69/NFAT reporter activation, cited to a specific paper or IEDB record). No entry is a prediction or model output — this is a literature/database curation exercise only.

## 2. Inclusion / exclusion criteria actually applied

**Cognate:** included only if (a) a solved PDB ternary structure (TCRα+TCRβ+MHC-IIα+MHC-IIβ+peptide, all 5 chains actually present — verified by fetching the RCSB entry/FASTA myself, not assumed from a title) was confirmed, OR (b) a specific paper/IEDB record reported a positive functional assay (proliferation, cytokine release, tetramer/SPR binding, CD69/NFAT reporter) for a *named* TCR against a *named* peptide-MHC-II. CDR3 sequences for every `solved_structure` row were parsed by me directly from the actual RCSB-deposited chain FASTA using the standard IMGT convention (from the conserved Cys to the conserved Phe/Trp of the **F**G-X-G/**W**G-X-G junction motif, inclusive), not copied from a paper abstract or recalled from memory.

**Non-cognate:** included only if the **same** TCR (same CDR3a/CDR3b), on the **same** restricting MHC allele, was explicitly tested against a **different** peptide and reported as non-stimulatory/non-binding in a cited source (SPR no binding, tetramer-negative, "no response"/"did not recognize" in a functional assay). I did **not** create a non-cognate row for "an untested peptide" or infer non-cognate status from "this wasn't the peptide the structure was solved with." Several strong *candidate* non-cognate pairs were found in the literature but excluded because the response was **reduced, not absent** (e.g., TCR T1D3 responds ~1000× more weakly, but not zero, to the wild-type insulin peptide vs. its mimotope) — a quantitative-but-nonzero difference is not the same as a documented negative, so it was left out rather than stretched into a non-cognate claim.

**Cross-reactive pairs are NOT the same as non-cognate pairs**, and this distinction mattered a lot in practice: several papers in this space (the GAD65/flu seed clones, the DQ2.5-gliadin XPA5 clone, the MS Ob.1A12/2D2 mouse clones, the HLA-DQ1 Hy.1B11 clone) report a **second peptide that also activates** the same TCR (molecular mimicry / degenerate recognition) — those second peptides were coded as **additional cognate rows**, not non-cognate rows, per the task's explicit definition.

**Provenance discipline:** every PDB ID was confirmed to exist and to have the claimed 5-chain composition via a live RCSB fetch in this session (not assumed from a paper title). Every non-cognate claim traces to a specific assay/PMID/IEDB record, cited in `source_id`/`source_url`. Where a sub-agent's extracted fact could not be independently reproduced, I flagged and corrected it rather than trusting it silently (see §5).

## 3. Class balance and allele skew

**Class balance:** 36 cognate / 5 non-cognate (88% / 12%) — heavily skewed toward cognate, because documented experimental *negatives* are far rarer in the published literature than documented positives (papers report what worked; explicit "we tested peptide X and got no response" statements are comparatively uncommon and often buried in supplementary text that isn't fully accessible via automated fetching).

**MHC allele skew — confirmed, and substantial:**

| MHC β allele | rows | % of 41 |
|---|---|---|
| HLA-DRB1\*04:01 | 12 | 29% |
| HLA-DQB1\*02:01 (DQ2.5) | 12 | 29% |
| HLA-DQB1\*03:02 (DQ8) | 6 | 15% |
| mouse I-Ag7 / I-Ak / I-Au / I-Ab (combined) | 7 | 17% |
| HLA-DRB1\*01:01 | 2 | 5% |
| HLA-DRB1\*15:01, HLA-DQB1\*05:02 | 1 each | 2% each |

**DRB1\*04:01 and DQ2.5 together account for 24/41 rows (59%) of the whole benchmark.** This mirrors the seed data (which is 100% DRB1\*04:01) and is not an artifact of my search strategy alone — it reflects a genuine, well-documented skew in the *solved-structure and IEDB literature itself*: DRB1\*04:01 (autoimmune/RA/T1D disease-association structural biology, Rossjohn/Reid/Kappler/Mariuzza labs) and DQ2.5 (celiac disease structural biology, the single most structurally over-represented class-II system in the PDB) dominate the field of solved class-II ternary TCR-pMHC structures. If this benchmark is used downstream, this skew should be corrected for (e.g., by down-weighting or by deliberately seeking out additional DP- and rarer-DQ-restricted entries) rather than treated as representative of class-II restriction generally.

**Species:** 32 rows are human MHC (with a few using mouse TCRs raised in HLA-DR4-transgenic mice, still human MHC-restricted); 9 rows are fully mouse (I-Ag7 ×3, I-Ak ×2/×2, I-Au ×1/×1, I-Ab ×1).

## 4. Gaps — what I looked for but could NOT verify

- **Wucherpfennig & Strominger 1995 (Cell, PMID 7534214)**, the foundational MBP/EBV molecular-mimicry paper — abstract says "seven viral and one bacterial peptide activated" MBP-specific clones "and 129 peptides were tested" on 7 clones, implying many explicit negatives exist, but the full peptide-by-peptide activation table was not retrievable via WebFetch (paywalled full text) — not included rather than guessed.
- **DQ2.5-glia-α1a/ω1 TRAV9-2+/TRBV7-3+ non-cognate pair** (Ciacchi 2022 JBC) — confirmed to exist in principle (this is exactly the W316/W321 pair I *did* include), so no gap there, but a companion "PB"-family citrullinated-Tenascin-C **single-citrullination-register** negative (P-1-only or P2-only citrulline, mentioned in Dao et al. 2025 as also reduced/negative) was not fully extracted with sequence-level detail and was left out.
- **Biased TCR gene usage in citrullinated Tenascin-C T cells** (PMC8720095, Sci Rep 2021) — confirms a same-TCR-tested-against-a-different-citrullinated-epitope negative ("cit-TNC22 TCR … not … cit-TNC17") exists, but the actual CDR3 sequences are in a supplementary Excel/table not retrievable via WebFetch — excluded rather than fabricated.
- **8F10 TCR (I-Ag7, insulin mimotope register-3B)** — secondary sources state it "responded ONLY to the register-3B tetramer," strongly implying an explicit negative against another tested register/native peptide, but the exact negative peptide sequence could not be extracted from accessible text in this session — flagged as a likely-real gap rather than invented.
- **TRAV/TRAJ/TRBV/TRBJ gene names** for roughly half of the `solved_structure` rows (all of the DQ2.5-gliadin Petersen-2014 clones JR5.1/D2/S16/S2, the DQ8 clones SP3.4/T1D3/A2.13, the mouse D10/172.10/4.1/8F10 clones, and the citrullinated-antigen clones A07/A03/RA2.7) could **not** be confirmed from an accessible primary-source table (paywalled full text, or supplementary tables not retrievable via WebFetch) even though the **CDR3 amino-acid sequences themselves are independently verified** (parsed directly from the deposited PDB FASTA). These are marked `NOT FOUND` in the TRAV/TRAJ/TRBV/TRBJ columns rather than guessed — usable for CDR3-only benchmarking, incomplete for full V/J-gene-level analysis.
- **MHC_alpha_seq / MHC_beta_seq** for the 5 IEDB-only rows (Ob.1A12, Hy.1B11, 2D2, BDC-2.5, TAZ10) are `NOT FOUND` — I have the allele *names* from IEDB but did not fetch the actual MHC chain sequences for these (no PDB structure exists for these particular clone/peptide pairs; would require an IMGT/HLA or UniProt lookup not done in this session).
- **PDB 8VCY and 6XCP** (A2.13 + HIP3 / HIP1 structures) — chain composition/peptide identity relayed from the primary paper's structure list by a sub-agent and from a citation respectively, not independently re-fetched/re-parsed by me from RCSB directly (unlike every other `solved_structure` row in this file, which I did verify myself) — flagged in the relevant `confidence_note` fields for follow-up.

## 5. A caught error — why every fact here was independently re-checked

One sub-agent's initial extraction of CDR3 sequences for PDB 8TRQ (TCR "A07") from a secondary source did **not** match what I obtained by parsing the actual deposited PDB FASTA myself. I did not use the sub-agent's unverified figures; I substituted my own directly-parsed values and left an explicit note on that row. Given this, I independently re-derived CDR3 sequences from RCSB FASTA myself (not from any sub-agent's summary) for every `solved_structure` row in the final CSV except 8VCY/6XCP (flagged above), and cross-validated two IEDB-only records (S2 and FS18) against their own solved structures, both of which matched exactly — giving reasonable but not absolute confidence in the remaining IEDB-only rows, which are flagged as such.

## 6. Seed data (`TCRmodel2_submission_sheet.csv`) — important reclassification

The three seed clones (Clone_81, Clone_566, P196-1) were traced to a real, confirmed primary publication: **Cerosaletti K et al., "Germline-like TCR-α chains shared between autoreactive T cells in blood and pancreas," Nature Communications 2024;15:4971 (PMID 38871688, DOI 10.1038/s41467-024-48833-w)**, plus an earlier related paper, Linsley PS et al., JCI Insight 2021;6(22):e151349 (PMID 34806648). The primary paper's own TCR-gene table matches the seed file's V/J genes and CDR3 sequences exactly.

**However**, the primary source states that **both** GAD65 377-396 **and** influenza MP54 97-116 peptides activated Clone_81, Clone_566, and P196-1 **identically** — i.e., these are multi-specific/cross-reactive TCRs (a genuine molecular-mimicry finding, which is directly on-topic for this project), not TCRs with one cognate and one non-cognate peptide. **All six seed-derived rows are therefore coded `cognate` in the final CSV; none is coded `non_cognate`**, contradicting the seed file's original cognate/non-cognate labels for the "opposite" peptide of each clone. This is flagged in each row's `confidence_note`. Also note: the seed file names the influenza peptide "MP74_97-116" while the primary publication names the same peptide "MP54 97-116" — likely a naming discrepancy between the local project file and the paper; the sequence (`VKLYRKLKREITFHGAKEIS`) was carried over unchanged from the seed file since I could not independently confirm which numbering is correct.

## 7. Files

- `tcr_pmhc_benchmark.csv` — 41 rows, columns: `pair_id, class, tcr_source_id, TRAV, TRAJ, CDR3a, TRBV, TRBJ, CDR3b, peptide_name, peptide_seq, MHC_class, MHC_alpha_allele, MHC_beta_allele, MHC_alpha_seq, MHC_beta_seq, evidence_type, source_id, source_url, confidence_note`.
- This notes file.
