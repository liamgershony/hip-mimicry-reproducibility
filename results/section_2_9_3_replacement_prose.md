# Phase 1.5 — draft replacement prose for §2.9.3

## Current text (as in main (4).pdf)

> **2.9.3. Composite prioritization**
> Each metric was expressed as percentage change relative to the native HIP and
> signed so that positive values indicated improvement:
> %Δ = (x_candidate − x_HIP) / x_HIP × 100.  (3)
> Within-system sample standard deviations were calculated as
> s_Δ = sqrt( Σ(Δᵢ − Δ̄_cohort)² / (n−1) ).  (4)
> where Δᵢ is the percentage change for candidate i and Δ̄_cohort is the mean
> percentage change across the system's candidate cohort... Each percentage
> change was divided by its within-system standard deviation, and the three
> standardized values were averaged with equal weights.

**The problem:** "the system's candidate cohort" is never defined, and the
natural reading — all Tier 1 candidates screened in that system (n=28 for
HIP6/A2.11/DQ8, n=20 for HIP11/8.E3/DQ8-trans) — is not what was computed.
Verified in Phase 1 (`scripts/build_v3_figures.py` diff and manual
recomputation, `results/` not yet written to a dedicated script but confirmed
by hand): the SDs in Eq. 4, and every composite score in Table 10, reproduce
exactly **only** when the cohort is the small set of candidates that already
passed the native-HIP CI gate on the primary criterion — n=4 for HIP6/A2.11/DQ8,
n=6 for HIP11/8.E3/DQ8-trans — not the full Tier 1 pool for that system.

## Suggested replacement

> **2.9.3. Composite prioritization**
> Each metric was expressed as percentage change relative to the native HIP and
> signed so that positive values indicated improvement:
> %Δ = (x_candidate − x_HIP) / x_HIP × 100.  (3)
> **The standardization cohort for Equation 4 was the set of candidates already
> meeting the native-HIP criterion in that system (Section 2.9.2) — four
> candidates for HIP6/A2.11/DQ8 and six for HIP11/8.E3/DQ8-trans — not the full
> Tier 1 candidate pool screened in that system. Within-system sample standard
> deviations were therefore calculated over this already-filtered cohort:**
> s_Δ = sqrt( Σ(Δᵢ − Δ̄_cohort)² / (n−1) ),  (4)
> where Δᵢ is the percentage change for candidate i and Δ̄_cohort is the mean
> percentage change across **this filtered cohort, not the full screened
> population**. Each percentage change was divided by its within-system
> standard deviation, and the three standardized values were averaged with
> equal weights. **Because the denominator is estimated from the same small,
> already-selected group being ranked, the composite score describes relative
> position within the winners rather than deviation from the broader candidate
> population, and its magnitude is not comparable to a standardization computed
> against the full Tier 1 pool.**

## Cross-reference

§3.15 already states: "Composite ranks depend on cohort membership. The
standardization value is calculated from the candidates being ranked. Adding or
removing a candidate can therefore change their order." This replacement moves
that existing admission — which currently only says the SD is "calculated from
the candidates being ranked" without stating what that cohort actually *is* in
absolute terms (n=4/n=6 vs. n=28/n=20) — into the Methods section where a
reader building on this work would need it before, not after, encountering
Table 10.
