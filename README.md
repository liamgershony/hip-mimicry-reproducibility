# Reproducibility archive — HIP structural mimicry screen

Data and analysis code for:

> **A two-tiered in silico pipeline prioritizes candidate viral and bacterial mimics
> of hybrid insulin peptide autoantigens in type 1 diabetes**
> Liam Gershony, James M. Heather

Every table, figure, and statistic reported in the manuscript regenerates from the
files here with one command.

```bash
./run_all.sh
```

Requires `python3` with `openpyxl`, `scipy`, `numpy`. Step 7 additionally needs
[`tectonic`](https://tectonic-typesetting.github.io/) to compile the TikZ figure; it
is skipped automatically if not installed.

---

## Provenance — read this first

Not all values in this archive have the same evidentiary status, and the difference
matters. Every per-seed value in `data/per_seed_metrics.csv` carries a `provenance`
tag:

| Tag | Meaning | Count |
|---|---|---|
| `RAW` | Traceable to an actual model-output file included here | 25 |
| `WORKBOOK` | Exists only as a curated value in the results workbook; no upstream per-job model output was retained | 1,220 |
| `MISSING` | Not recorded anywhere | 610 |

Concretely:

- The **five reference-complex AlphaFold 3 jobs** are genuine raw model output —
  CIF structures, per-seed `summary_confidences` JSON, and the `job_request.json`
  recording the exact submitted sequences. These are in
  `alphafold3_reference_complexes/`.
- The **55 Tier 1 candidates, 3 negative controls, and 3 native-HIP re-runs** exist
  only as curated spreadsheet records. The original per-job TCRModel2 and PRODIGY
  output files were not retained and are therefore not in this archive. Reanalysis
  of the curated records cannot recover them.
- **Overall model pLDDT, pTM, and model-confidence** were never recorded for any of
  the 61 peptides. They are represented as explicit `MISSING` rows rather than
  omitted, so the gap is visible in the data itself.

This is stated plainly because the alternative — implying a completeness the record
does not support — would be worse. `results/` contains the full audit trail.

---

## Layout

```
data/
  Master_Results_Sheet.xlsx          source workbook (the upstream record)
  per_seed_metrics.csv               tidy: one row per system/peptide/seed/metric,
                                     with provenance + originating workbook cell
  per_peptide_summary.csv            per-peptide stored RMSD, p-value, CI, SD
  reference_complex_af3_metrics.csv  the 25 RAW ipTM/pTM values
  tcr_pmhc_benchmark.csv             36 cognate + 5 documented non-cognate pairs
  tcr_pmhc_benchmark_decoys.csv      100 shuffled decoys (assumed, not verified)
  tcrmodel2_submission_sheet.csv     202 planned jobs, inputs + empty metric columns
  prodigy_worklist.csv               1,010 planned PRODIGY runs (202 jobs x 5 models)
  *_notes.md                         curation criteria, class balance, allele skew

scripts/                             all analysis code (paths are repo-relative)
results/                             every regenerated output
corrections/                         flagged data issues, uncorrected by design
alphafold3_reference_complexes/      raw AlphaFold 3 output, 5 systems x 5 models
```

## What the scripts reproduce

| Script | Reproduces |
|---|---|
| `build_csv_spine.py` | The tidy CSV spine + provenance tags from the workbook |
| `recompute_welch.py` | The 48/55 Welch's t-test count |
| `recompute_ci_gate.py` | 12 candidates passing the four-metric native-HIP gate; 10 after the human-association filter |
| `reanalysis_no_new_data.py` | Mann–Whitney U, Benjamini–Hochberg FDR, Cliff's delta with bootstrap intervals |
| `kd_dg_consistency.py` | Kd/ΔG internal consistency across all 305 seed values |
| `derive_tolerance_band.py` | The pre-declared reproduction tolerance band |
| `rebuild_figure3.py` | Figure 3, rebuilt from verified data |
| `diff_figure3_transcription.py` | Audit of an earlier screenshot-derived figure against the workbook |
| `build_submission_sheet.py`, `build_prodigy_worklist.py`, `build_shuffled_decoys.py` | Planned-work sheets (not yet executed) |

## Known data issues, deliberately not silently corrected

- `corrections/kd_dg_outliers_FOR_REVIEW.csv` — one negative-control cell records
  Kd = 1.4 × 10⁻¹⁰⁷ against a ΔG of −9.4 kcal/mol, which implies ~1.3 × 10⁻⁷. It is
  an isolated data-entry error, not systemic: across all 305 seed values the median
  |ΔG − RT·ln(Kd)| is 0.024 kcal/mol. It does not propagate into any reported table.
- `results/figure3_transcription_diff.csv` — 15 of 275 ipTM values in an earlier
  screenshot-derived version of Figure 3 disagreed with the workbook. That figure was
  retracted and rebuilt from source; the diff is kept as an audit record.

## Planned but not executed

`tcrmodel2_submission_sheet.csv` and `prodigy_worklist.csv` specify a full
regeneration and a cognate/non-cognate benchmark that were **not run**. They are
included because they document exactly what would be required, including the PRODIGY
chain selection (`D,E` TCR vs `A,B,C` pMHC) derived from the TCRmodel2 source. No
result in the manuscript depends on them.

## Citation

Please cite the manuscript. If you use the benchmark curation, note that the 100
shuffled decoys are *assumed* negatives based on TCR specificity, not experimentally
verified — see `data/tcr_pmhc_benchmark_decoys_notes.md`.
