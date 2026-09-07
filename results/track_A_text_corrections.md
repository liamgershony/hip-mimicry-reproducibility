# Track A — text corrections (drafts for Prism / manuscript editing)

## A1. §3.7 replacement prose

**Current text (main (4).pdf, §3.7):**
> Ten of the 48 remaining candidates met or exceeded all four thresholds: four
> in HIP6/A2.11/DQ8, none in HIP11/8.E3/DQ8, and six in HIP11/8.E3/DQ8-trans.
> Their organism-level distribution is summarized in Table 7.

**Problem:** this states the CI gate alone produced 10. It produced 12
(`Reproducibility_Pipeline/results/ci_gate_recompute.csv`); 10 remained only
after the separately-documented human-association filter (§2.4) was applied.

**Replacement:**
> Twelve of the 48 remaining candidates met or exceeded all four native-HIP
> thresholds on the structural and energetic criteria (Section 2.9.2). Two of
> these twelve — both animal-associated Rotavirus A strains (Section 3.6) —
> were then excluded under the human-association nomination criterion
> (Section 2.4), which is applied at the nomination stage rather than as part
> of the four-metric filter itself. Ten candidates remained after both steps:
> four in HIP6/A2.11/DQ8, none in HIP11/8.E3/DQ8, and six in
> HIP11/8.E3/DQ8-trans. Their organism-level distribution is summarized in
> Table 7.

**Figure 1 flowchart change:** insert a new box between the "Native-HIP 95% CI
gate" node and the three-lane 4/0/6 split:

```
Native-HIP 95% CI gate
        |
        v
  "12 pass all four metrics"          <- NEW intermediate node
        |
        v
Human-association filter              <- NEW node
(excludes 2 animal-associated
 Rotavirus A candidates, Sec 2.4)
        |
        v
   [4]      [0]      [6]   <- existing 3-lane split, now correctly downstream
        \    |    /
         v   v   v
     "10 candidates"
```
Implemented in `Figures/latex_source/fig_pipeline.tex` (see diff below); Prism
should re-place this as the new Figure 1 once approved.

---

## A2. §2.6 replacement prose

**Current text (main (4).pdf, §2.6):**
> Recorded metrics included pLDDT, pTM, ipTM, TCR-pMHC ipTM, interface-pLDDT,
> and model confidence.

**Problem, verified against the actual workbook column headers** (32 columns,
enumerated in `scripts/build_csv_spine.py`): only **TCR-pMHC ipTM** and
**interface-pLDDT** exist as recorded per-seed metrics for the 55 Tier 1
candidates. Overall model pLDDT and pTM do not exist as columns at all.
**Additional finding while drafting this correction**: a separate local file
(`~/Downloads/HIP11_8E3_TRB_ranking.csv`) shows "model confidence" (as
`Mean_ModelConf`) *was* recorded at an earlier pipeline stage — TCR
beta-chain selection for clone 8.E3 — but that column does not exist in
`Master Results Sheet.xlsx` for the 55 final candidates either. So "model
confidence" is unsupported for the reported results in the same way pLDDT and
pTM are, even though the concept was used upstream. This should be added to
Track A3's provenance table as a third MISSING metric, with this same caveat
(recorded at TCR-selection stage, not retained for final candidate results).

**Replacement:**
> Recorded metrics were TCR–pMHC ipTM and interface-pLDDT, per seed. Overall
> model pLDDT and pTM, and a separate per-candidate model-confidence score
> used during earlier TCR beta-chain selection, were not retained in the
> results workbook for the 55 Tier 1 candidates (Supplementary Table SX).

---

## A3. Supplement provenance section

**Draft, written to be revised upward if Track C succeeds:**

> ### Supplementary Table SX. Data provenance
>
> | Provenance level | Definition | Count (per-seed metric values) |
> |---|---|---|
> | RAW | Traceable to an actual model-output file on disk | 25 (the 5 reference-complex systems' AlphaFold 3 job outputs only) |
> | WORKBOOK | Exists only as a value in the results workbook; no upstream model-output file has been located | 1,220 (TCR–pMHC ipTM, interface-pLDDT, ΔG, and Kd for the 55 Tier 1 candidates, 3 negative controls, and 3 TCRModel2-based native-HIP re-runs) |
> | MISSING | Not present anywhere | 610 confirmed (overall pLDDT and pTM for all 61 peptides); model-confidence scores for the same 61 peptides are also absent from the results workbook, though the concept was recorded at an earlier, separate TCR-selection stage not part of the final candidate results |
>
> Of the two datasets this study depends on, only the five reference-complex
> models (Table 2, Figure 2) are backed by raw, machine-readable model output.
> The 55-candidate Tier 1 dataset that Tables 8, 9, and 10 and Figure 3 are
> built from exists only as curated, two-decimal-precision values in a results
> spreadsheet; no per-seed job files, model version strings, random seeds, or
> command lines for these runs have been located. This dataset should be
> described as spreadsheet-level, not as a reproducibility package, until raw
> per-seed outputs are regenerated and archived (Section 2.10).
>
> *This statement will be revised if the pilot regeneration described in
> Section 2.10 succeeds and a full rerun is subsequently completed and
> archived — see the corresponding update to the Data Availability statement.*

---

## Status
All three drafts above are ready for Prism to apply. None have been applied to
the manuscript yet — reporting per the gate instruction.
