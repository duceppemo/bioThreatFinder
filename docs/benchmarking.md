# Benchmarking

Every threshold in `packs/*/pack.yaml` started as a reasoned guess. This harness
is what replaces guesses with measurements.

## Read this first: what simulation cannot tell you

Simulated reads flatter every tool. Error profiles are idealised, there is no
real reagent contamination, no index-hopping chemistry, no degraded or inhibited
clinical material, no PCR duplication artefacts, and — most importantly — **no
genuinely novel agent**, because anything in a database to simulate from is by
definition not novel.

So these scenarios establish **relative** performance and catch regressions.
They do **not** establish clinical performance. That needs, in rough order of
increasing value:

1. Mock communities with certified composition (e.g. ZymoBIOMICS)
2. Contrived specimens — spiked matrix at known titres, processed through the
   real wet-lab protocol so extraction and host depletion losses are included
3. Proficiency-test panels
4. Characterised real specimens with an orthogonal reference method

Until at least stage 2 is done, no accuracy claim about patient specimens is
assertable. The harness will happily print a beautiful LOD that the wet lab
cannot reproduce.

**Two limits no amount of benchmarking removes:**

- **Software cannot beat its input.** For low-biomass specimens the binding
  constraint is host fraction, extraction efficiency and sequencing depth, not
  the classifier. That is why `lod_lowbiomass_ypestis` runs at 99.9% host — the
  number from that scenario, not the 95% one, is the honest practical LOD.
- **Homology-based detection has a ceiling.** A wholly novel or engineered agent
  with no detectable homology to anything known is invisible to every tier here,
  including HMMs. `divergent_*` scenarios probe the edge of that envelope using
  untargeted family members as stand-ins; they do not abolish it.

## Design: scenarios built to break the tool

A benchmark that spikes targets into a clean background measures almost nothing —
any tool looks good at 100 RPM of pure *B. anthracis*. The 25 scenarios in
`benchmarks/scenarios/human_core.yaml` are weighted towards failure modes:

| Kind | n | What it measures |
|---|---|---|
| `lod_series` | 6 | Sensitivity and LOD across genome sizes (30 kb virus → 23 Mb parasite), plus a 99.9%-host case |
| `specificity` | 9 | The near neighbours that actually cause false positives, with the target genuinely **absent** |
| `marker_ablation` | 3 | Target chromosome **without** its virulence loci |
| `vaccine` | 4 | BCG, vaccinia, Sabin, RB51 — must read as stand-down, not disease |
| `negative` | 1 | Host + commensals only; sets the run-level false-positive rate |
| `contamination` | 2 | Index hopping into a negative, and a batch with no control at all |

The `exclude_loci` mechanism is what makes the marker scenarios possible: it
simulates *B. anthracis* chromosome with pXO1/pXO2 **removed**, so alignment
evidence is overwhelming and the correct answer is still "not confirmed".

23 expectations are machine-checkable, so `evaluate_benchmark.py` exits non-zero
on regression and can gate CI.

## Semantics, stated explicitly

Tiers map to four outcomes. These are defined in code, not left implicit:

| Outcome | Tiers |
|---|---|
| POSITIVE | `supported`, `confirmed_in_silico` |
| VACCINE | `vaccine_strain` |
| NEGATIVE | `screened_only`, `attributed_to_neighbour`, `marker_veto`, `suppressed_control`, `control_background`, or no finding |
| NOT_ASSESSED | `no_reference` |

A fifth case: an organism spiked as **background or as a near-neighbour confuser
that is not on the watchlist**. The pipeline is not looking for it, so scoring it
as a detection target would print "LOD not reached" for an organism no LOD exists
for. Pass `--watchlist` and these are marked `present_untargeted` and excluded
from sensitivity and LOD. Without `--watchlist` the evaluator warns.

**`NOT_ASSESSED` is excluded from every numerator and denominator.** If a missing
reference genome counted as a true negative, adding watchlist taxa without
building their references would appear to *improve* specificity. The
self-test asserts this exclusion.

## Metrics

- **Sensitivity** per taxon × platform × spike level, across replicates
- **LOD** — lowest level where *every* replicate is POSITIVE (`lod_all_rpm`), plus
  the 50% point. Flagged `coarse` below 5 replicates, because it is
- **Specificity** — POSITIVE and `confirmed_in_silico` rates on taxa the truth
  says are absent. Weighted most heavily in the summary: a false confirmation on
  a high-consequence organism triggers a real public-health response
- **Tier PPV** — the calibration number. If `confirmed_in_silico` does not show
  markedly higher PPV than `supported`, the tier structure is not earning its
  complexity and should be simplified
- **Vaccine accuracy** — fraction correctly tiered rather than reported as disease
- **Ablation** — the same data scored with and without marker gating, quantifying
  what the marker tier actually contributes

## Running it

```bash
# 1. Generate datasets (needs InSilicoSeq/ART and Badread/NanoSim)
python3 scripts/simulate_benchmark.py \
    --scenarios benchmarks/scenarios/human_core.yaml \
    --batch v1 --outdir benchmarks/runs --seed 42

# 2. Run the pipeline over the generated samplesheet
nextflow run . -profile conda,workstation \
    --input benchmarks/runs/v1/samplesheet.csv \
    --pack human --outdir benchmarks/runs/v1/results

# 2b. Optional ablation arm: identical data, marker gating disabled
nextflow run . -profile conda,workstation \
    --input benchmarks/runs/v1/samplesheet.csv \
    --pack human --skip_markers \
    --outdir benchmarks/runs/v1/results_nomarkers

# 3. Measure
python3 bin/evaluate_benchmark.py \
    --scenarios benchmarks/scenarios/human_core.yaml \
    --manifest benchmarks/runs/v1/manifest.tsv \
    --watchlist packs/human/watchlist.resolved.tsv \
    --findings benchmarks/runs/v1/results/score_evidence/findings.tsv \
    --ablation-findings benchmarks/runs/v1/results_nomarkers/score_evidence/findings.tsv \
    --out-dir benchmarks/runs/v1/evaluation
```

Outputs: `summary.md`, plus `sensitivity.tsv`, `lod.tsv`, `specificity.tsv`,
`tier_ppv.tsv`, `vaccine.tsv`, `expectations.tsv`, `ablation.tsv` and the raw
`observations.tsv`.

## Self-tests

Both stdlib-only and runnable with no data:

```bash
python3 bin/btf_yaml.py --self-test           # config parser
python3 bin/marker_rules.py --self-test       # three-valued rule logic
python3 bin/evaluate_benchmark.py --self-test --scenarios x --manifest x --findings x --out-dir x
```

## A note on the background composition

`defaults.background` puts 2% commensal and environmental organisms in every
library. Without it, "negative" means host-only, which understates the
false-positive rate — real specimens are full of organisms genomically adjacent to
watchlist taxa.

Every background taxon is a near neighbour already in the pack's reference
manifest, so nothing extra needs fetching. Note *E. coli* in particular: commensal
*E. coli* **should** be detected and should **not** be confirmed, because
confirmation is gated on pathotype markers it does not carry. Having it in every
background is a free, continuous test of that rule.

Background taxa are deliberately absent from per-scenario `truth` blocks, so they
report as `unlisted` and stay out of the PPV denominator. `negative_host_only`
lists them explicitly, so that one scenario's false-positive rate is exact.

## Tool bake-off: comparing tools instead of arguing about them

Tool choice in this field is usually argued from papers. Papers benchmark on their
authors' data, with their authors' thresholds, usually for a community-profiling
question rather than an "is organism X present" question. The only answer that binds
for *this* pipeline, *this* watchlist and *this* specimen type is a measurement on
identical input.

`scripts/run_bakeoff.sh` runs the pipeline once per configuration over the **same**
simulated batch and compares them. The variant matrix is deliberately built from
**minimal pairs** — a variant that changes two things at once cannot attribute the
difference to either, which is the commonest way a tool comparison fails to answer
the question it was run for.

| Variant | Isolates |
|---|---|
| `kraken2_minimizer` | Baseline: Kraken2 `--report-minimizer-data` + Bracken + sylph |
| `plus_krakenuniq` | **Does HyperLogLog cardinality beat Kraken2 minimizer counts?** |
| `plus_metabuli` | Value of the DNA+amino-acid arm (watch the `divergent_*` scenarios) |
| `all_arms` | Upper bound on sensitivity *and* on cost |
| `no_sylph` | Is ANI/containment contributing anything the k-mer arms don't? |
| `no_markers` | What marker gating buys |

```bash
scripts/run_bakeoff.sh --batch v1 --list      # show the matrix
scripts/run_bakeoff.sh --batch v1 --dry-run   # show the commands
scripts/run_bakeoff.sh --batch v1
```

It is resume-friendly: a completed variant is skipped, and one failing variant is
recorded and omitted rather than abandoning the other five.

### Cost is half the comparison

The bake-off table reports CPU-minutes and peak RSS alongside accuracy, harvested
from each run's Nextflow trace. This is not decoration. KrakenUniq costs roughly
1400 GiB transient disk, a 250–450 GiB resident index, a **Jellyfish 1.x** build
dependency (this box has 2.3.1, which will fail) and 12+ hours to build — while
Kraken2's `--report-minimizer-data` already yields distinct minimizer counts and a
duplicity ratio from an index we build anyway. A configuration that gains a point of
sensitivity for that has to justify it, and a comparison that hides the cost side
invites exactly the wrong decision.

Results are ordered by **fewest false confirmations first**, then sensitivity: on a
biothreat watchlist a false confirmation triggers a real-world response, so it costs
more than a missed borderline lead.

### Read the per-taxon table before ranking anything

`bakeoff_per_taxon.tsv` gives detection rate, confirm rate and LOD per taxon *per
variant*. An aggregate sensitivity figure hides the case where one configuration
finds the divergent filovirus and another finds the extra *Plasmodium*. That is a
trade to decide deliberately, not a ranking to read off.

## Once results exist

Feed them back:

1. Set `min_breadth`, `min_unique_reads` and `min_occupancy` per taxon from the
   measured curves instead of by reasoning
2. Where a taxon's LOD is unacceptable, decide explicitly whether to trade
   specificity for it — and record the decision in the pack
3. Replace fixed thresholds with an empirical null from a background panel, and
   report calibrated FDR (roadmap item 2)
4. Commit the results so regressions are visible in review
