# Track C1 — pricing, and Track C2 — pilot specification

## C1. Real-dollar cost estimates

### The most important finding: TCRModel2 and PRODIGY both have official free public web servers

- **TCRModel2**: `https://tcrmodel.ibbr.umd.edu/` (Pierce Lab, University of Maryland IBBR)
  — the actual web server behind the manuscript's own citation (Yin et al. 2023, NAR).
  Confirmed live via WebFetch. No published version number or queue-time SLA on the page
  itself. Free. This is almost certainly the closest available version-match to whatever
  the original 55-candidate runs used, since nothing in this workspace suggests a private
  GPU setup was ever used — the same pattern already holds for the reference complexes,
  which were run through the official AlphaFold Server (`alphafoldserver.com`), not a
  private install.
- **PRODIGY**: `https://wenmr.science.uu.nl/prodigy/` and `https://bianca.science.uu.nl/prodigy/`
  (Bonvin Lab, Utrecht University). Free, official, matches the manuscript's citation.
- **Implication for the pilot (C2)**: the pilot should default to submitting through
  these two free official servers rather than renting a GPU, both because it costs
  nothing and because it's the best available version-match. Cloud GPU rental is the
  fallback if the free servers can't handle the pilot's volume responsively, and is the
  likely necessity for the full 61+benchmark batch (below), where submitting 100+ jobs to
  a shared academic server risks being both slow and inconsiderate to other users.

### Cloud GPU options, if self-hosting becomes necessary (verified via live search, sources below)

| Provider | GPU | On-demand rate | Notes |
|---|---|---|---|
| RunPod (Community Cloud) | A100 | **$1.39/hr** | Cheapest verified A100; billed per-second |
| RunPod (Secure Cloud) | A100 PCIe | $1.39–1.49/hr | Vetted datacenters, ~2x community for RTX 4090 |
| Lambda Labs | A10 | **$1.29/hr** | On-demand, no long-term commitment |
| Lambda Labs | A100 40GB | $1.99/hr | |
| AWS EC2 | g5.xlarge (A10G 24GB) | **$1.01/hr** (us-east-1) | Standard on-demand |
| AWS EC2 | g6.xlarge (L4) | $0.80/hr | Newer, better inference $/perf per some sources |
| Google Colab | T4 (free tier) | **$0** | 16GB VRAM, ~90min idle timeout / 12hr session cap — usable for the pilot's individual jobs, not a reliable batch-scale option |

None of these require heavy setup for a ColabFold-style pipeline specifically — the one
real archived job I found on disk (`Shorter_HIP_test_Using_TCRMODEL2_SEQ_cf30e.result.zip`)
is a standard `alphafold2_multimer_v3` ColabFold run, which is what all of these providers
support with a stock ColabFold install (no custom infrastructure needed).

### Real-world throughput baseline (not invented — from that one archived log)

That job (a 404-residue TCR–pMHC-II complex, comparable in scale to this paper's systems)
took **~18 minutes of GPU compute for 5 models** (154–281 s/model) plus ~2 min queue/MSA
time. This is the only empirical timing data available; I'm using it as the working
estimate, not asserting it represents every system's runtime.

### Dollar estimates

**Pilot (5 systems × 5 seeds, Track C2 below), if self-hosted instead of using free servers:**
~5 × 18 min ≈ 90 min compute → **$1.50–2.10** on any of the above providers. Trivial
either way — the pilot's cost is not a real constraint regardless of path chosen.

**Full-scale regeneration** — 61 systems (55 candidates + 3 controls + 3 native re-runs)
plus Track B's benchmark (~40–50 additional systems once curated; using ~45 as a
placeholder pending Track B's actual count):

| Seeds | Total systems | Sequential GPU time | RunPod Community A100 | Lambda A10 | AWS g5.xlarge |
|---|---|---|---|---|---|
| 5 | 61 + ~45 ≈ 106 | ~31.8 hr | ~$44 | ~$41 | ~$32 |
| 20 | 61 + ~45 ≈ 106 | ~127 hr (~5.3 days) | ~$177 | ~$164 | ~$128 |

**Caveats, stated plainly:**
- These extrapolate linearly from ONE observed job. Benchmark-set systems (TCR + pMHC-II
  only, no extra chains) are likely smaller/faster than the paper's full complexes, so this
  is probably a mild overestimate for that portion — not corrected for, since I have no
  empirical data for that system size specifically.
- Sequential time can be cut by running multiple rented instances in parallel; dollar cost
  does not change by doing so (same total GPU-hours), only wall-clock time does.
- **Dollar cost is not the binding constraint here** — even the 20-seed full-scale run is
  under $200 on any of these providers. The real gating factors are (a) wall-clock time,
  (b) whether a self-hosted rerun's TCRModel2/ColabFold version matches whatever the
  original runs used (unknown — Phase 0's core finding), and (c) research-server etiquette
  if using the free official servers at scale.

**Sources:** [Lambda Labs pricing – ComputePrices.com](https://computeprices.com/providers/lambda), [Lambda Labs – SynpixCloud](https://www.synpixcloud.com/blog/lambda-labs-gpu-pricing-2026), [RunPod pricing – Northflank blog](https://northflank.com/blog/runpod-gpu-pricing), [RunPod RTX 4090 page](https://www.runpod.io/gpu-models/rtx-4090), [AWS g5.xlarge – Vantage](https://instances.vantage.sh/aws/ec2/g5.xlarge), [AWS G5/G6 pricing guide – Wring](https://wring.co/blog/aws-gpu-instance-pricing-guide), [ColabFold Colab hardware discussion](https://saturncloud.io/blog/whats-the-hardware-spec-for-google-colaboratory/), [TCRmodel2 web server](https://tcrmodel.ibbr.umd.edu/), [TCRmodel2 paper, NAR](https://academic.oup.com/nar/article/51/W1/W569/7151345), [PRODIGY web server, Bonvin Lab](https://wenmr.science.uu.nl/prodigy/manual).

---

## C2. Pilot specification

### The 5 pilot systems

Per your spec — 3 candidates spanning the score range, plus 1 native HIP, plus 1 negative
control, all at 5 seeds, matching the paper's configuration as closely as determinable
(same TCR/MHC/peptide chain composition per `per_seed_metrics.csv`; TCRModel2's own default
settings otherwise, since no non-default parameters were ever recorded — Phase 0 finding):

| Role | System | Peptide | Sequence | Published mean TCR–pMHC ipTM |
|---|---|---|---|---|
| High-scoring nominated candidate | HIP11/8.E3/DQ8-trans | RRNVATLQAENVTG (SARS-CoV-2, lead candidate) | — | 0.842 |
| Mid-scoring candidate | HIP11/8.E3/DQ8-trans | REDTVSVKSEPVSE (HHV-5/CMV) | — | 0.696 (passed Welch but failed CI gate — genuinely mid-range, not cherry-picked) |
| Gate-failed candidate | HIP6/A2.11/DQ8 | GVEDVYLAGALEAQ (HHV-8) | — | 0.676 (failed the Welch test itself, p=0.4072) |
| Native HIP | HIP11/8.E3/DQ8-trans | SLQPLALEAEDLQV (native HIP11) | — | 0.822 |
| Negative control | HIP11/8.E3/DQ8-trans | KDVDAAVDAEVVQF (EBV AG876) | — | 0.258 |

All five drawn from the same system (HIP11/8.E3/DQ8-trans) except the gate-failed
candidate, so the pilot also directly tests whether a fresh 5-seed native/control pair
reproduces close to the *same* native-HIP and control means already used throughout
Tables 8–10 — a second, independent check beyond the three individual candidates.

### Pre-declared tolerance band (derived from data, not chosen)

From `scripts/derive_tolerance_band.py`, run against all 61 peptides' actual per-seed SDs
in `per_seed_metrics.csv` (full derivation and caveats in `results/pilot_tolerance_band.md`):

- **TCR–pMHC ipTM: published mean ± 0.036** (median-case). For the two systems above with
  known higher seed variance (the mid and gate-failed candidates), the conservative
  (90th-percentile) band **± 0.37** applies instead — flagged now, before seeing any
  pilot data, so this isn't a post-hoc excuse.
- **interface-pLDDT: published mean ± 2.65** (median-case); conservative band ± 15.7 for
  the higher-variance systems.

This is declared **before** any pilot run, exactly as instructed. A pilot value landing
outside the applicable band is a disagreement, full stop — reported as such per your
explicit instruction on how to handle non-reproduction.

---

## C3. What the pilot must record (the template for a full rerun)

For each of the 5 systems, capture and archive:
- Full command line or exact web-server submission parameters (sequences, chain order,
  any non-default settings — there should be none, since the paper never states any beyond
  Amber relaxation)
- The random seed (or seed range) used for each of the 5 models
- TCRModel2 version / git commit or web-server date-stamp, and ColabFold/AlphaFold-Multimer
  backend version if disclosed
- Per-seed raw JSON output (or equivalent) with ipTM, interface-pLDDT, and any other
  reported confidence metrics
- Per-seed PDB structure files
- PRODIGY's raw output log (ΔG, Kd) per selected model, plus the exact PRODIGY version/date
- Wall-clock time per job (for validating or correcting the Track C1 cost model)

This becomes `Reproducibility_Pipeline/data/pilot/` with one subfolder per system, and a
`pilot_manifest.csv` indexing all of the above — the actual template a full 61+benchmark
regeneration would replicate.
