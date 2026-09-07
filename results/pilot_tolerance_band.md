# Pilot reproduction tolerance band (Track C2)

Derived from the 61 peptides' actual per-seed SDs in per_seed_metrics.csv, not chosen a priori.

- TCR_pMHC_ipTM: n=61 peptides, median per-peptide SD=0.0207, 90th-pct SD=0.2103 -> tolerance half-width median-case=+/-0.0364, conservative(90th-pct)-case=+/-0.3692
- interface_pLDDT: n=61 peptides, median per-peptide SD=1.5110, 90th-pct SD=8.9497 -> tolerance half-width median-case=+/-2.6529, conservative(90th-pct)-case=+/-15.7130

## Recommended bands for the pilot

- **TCR–pMHC ipTM: published mean ± 0.036** (median-case). If the pilot system is a known erratic one (e.g., a FAIL-status candidate with high seed variance), use the conservative band instead: ± 0.369.
- **interface-pLDDT: published mean ± 2.65** (median-case); conservative band ± 15.71.

## Caveats (stated, not hidden)

1. This assumes the rerun's own seed-to-seed variability is similar in magnitude to the published data's -- an assumption, not a fact, since we don't know if the original runs and a rerun will show comparable stochastic spread. The pilot's own 5 fresh seeds will let us check this directly (compare the rerun's own SD against the published SD for the same system once the pilot exists).
2. t_crit=2.776 is carried over from the paper's own CI convention (df=4, one 5-seed sample) for consistency; the technically correct df for a *difference* of two independent 5-seed means is smaller (Welch-Satterthwaite, typically ~6-8 when variances are similar), which would give a slightly larger t_crit and thus a slightly wider band than reported here. Reported bands are therefore mildly optimistic (narrower) rather than mildly conservative.
3. Per-peptide SD is highly right-skewed (mean >> median: 0.063 vs 0.021 for ipTM) -- driven by FAIL-status candidates with erratic, near-boundary seeds. The pilot spec should pick its 'mid' and 'fail' candidates aware that their own individual SD may run well above the population median, in which case the conservative (90th-percentile) band is the more honest one to apply to those specific systems.
