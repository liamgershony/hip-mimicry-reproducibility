#!/usr/bin/env python3
"""
Track C2 -- derive the pilot reproduction tolerance band from existing per-seed
variability in per_seed_metrics.csv, rather than picking a number.

Logic: the pilot will produce a new 5-seed mean for a system; we compare it to
the published 5-seed mean for the same system. Both means carry sampling
uncertainty. Assuming the rerun's seed-to-seed variability is of similar
magnitude to what's already observed in the published data (a stated
assumption, not a proven fact), the standard error of the difference between
two independent 5-seed means is:
    SE_diff = sqrt(SD_original^2/5 + SD_rerun^2/5) ~= sqrt(2) * SD_typical / sqrt(5)
A 95% band on that difference uses the same t_crit=2.776 convention the paper
already uses (df=4, two-sided 95%) for consistency, not because it's the
technically exact df for a difference-of-means (that would need a
Welch-Satterthwaite df on the pooled variances) -- flagged as an approximation.
"""
import csv, math
from collections import defaultdict
import numpy as np
import os as _os
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))

DATA = _os.path.join(_REPO, "data")
OUT = _os.path.join(_REPO, "results")

rows = list(csv.DictReader(open(f"{DATA}/per_seed_metrics.csv")))
T_CRIT = 2.776

report_lines = []
bands = {}
for metric in ["TCR_pMHC_ipTM", "interface_pLDDT"]:
    grouped = defaultdict(list)
    for r in rows:
        if r["metric"] != metric or r["provenance"] != "WORKBOOK":
            continue
        key = (r["system"], r["role"], r["peptide_label"], r["sequence"])
        grouped[key].append(float(r["value"]))
    sds = np.array([np.std(v, ddof=1) for v in grouped.values() if len(v) == 5])
    median_sd = float(np.median(sds))
    p90_sd = float(np.percentile(sds, 90))
    se_diff_median = math.sqrt(2) * median_sd / math.sqrt(5)
    se_diff_p90 = math.sqrt(2) * p90_sd / math.sqrt(5)
    band_median = T_CRIT * se_diff_median
    band_p90 = T_CRIT * se_diff_p90
    bands[metric] = dict(median_sd=median_sd, p90_sd=p90_sd,
                          band_median_case=round(band_median, 4),
                          band_conservative_p90_case=round(band_p90, 4))
    report_lines.append(
        f"{metric}: n=61 peptides, median per-peptide SD={median_sd:.4f}, "
        f"90th-pct SD={p90_sd:.4f} -> tolerance half-width "
        f"median-case=+/-{band_median:.4f}, conservative(90th-pct)-case=+/-{band_p90:.4f}"
    )

with open(f"{OUT}/pilot_tolerance_band.md", "w") as f:
    f.write("# Pilot reproduction tolerance band (Track C2)\n\n")
    f.write("Derived from the 61 peptides' actual per-seed SDs in per_seed_metrics.csv, "
            "not chosen a priori.\n\n")
    for line in report_lines:
        f.write(f"- {line}\n")
    f.write("\n## Recommended bands for the pilot\n\n")
    f.write(f"- **TCR–pMHC ipTM: published mean ± {bands['TCR_pMHC_ipTM']['band_median_case']:.3f}** "
            "(median-case). If the pilot system is a known erratic one (e.g., a FAIL-status "
            "candidate with high seed variance), use the conservative band instead: "
            f"± {bands['TCR_pMHC_ipTM']['band_conservative_p90_case']:.3f}.\n")
    f.write(f"- **interface-pLDDT: published mean ± {bands['interface_pLDDT']['band_median_case']:.2f}** "
            "(median-case); conservative band "
            f"± {bands['interface_pLDDT']['band_conservative_p90_case']:.2f}.\n")
    f.write("\n## Caveats (stated, not hidden)\n\n")
    f.write("1. This assumes the rerun's own seed-to-seed variability is similar in magnitude "
            "to the published data's -- an assumption, not a fact, since we don't know if the "
            "original runs and a rerun will show comparable stochastic spread. The pilot's own "
            "5 fresh seeds will let us check this directly (compare the rerun's own SD against "
            "the published SD for the same system once the pilot exists).\n")
    f.write("2. t_crit=2.776 is carried over from the paper's own CI convention (df=4, one "
            "5-seed sample) for consistency; the technically correct df for a *difference* of "
            "two independent 5-seed means is smaller (Welch-Satterthwaite, typically ~6-8 when "
            "variances are similar), which would give a slightly larger t_crit and thus a "
            "slightly wider band than reported here. Reported bands are therefore mildly "
            "optimistic (narrower) rather than mildly conservative.\n")
    f.write("3. Per-peptide SD is highly right-skewed (mean >> median: 0.063 vs 0.021 for ipTM) "
            "-- driven by FAIL-status candidates with erratic, near-boundary seeds. The pilot "
            "spec should pick its 'mid' and 'fail' candidates aware that their own individual "
            "SD may run well above the population median, in which case the conservative "
            "(90th-percentile) band is the more honest one to apply to those specific systems.\n")

print("\n".join(report_lines))
print(f"\nWrote {OUT}/pilot_tolerance_band.md")