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

## Once results exist

Feed them back:

1. Set `min_breadth`, `min_unique_reads` and `min_occupancy` per taxon from the
   measured curves instead of by reasoning
2. Where a taxon's LOD is unacceptable, decide explicitly whether to trade
   specificity for it — and record the decision in the pack
3. Replace fixed thresholds with an empirical null from a background panel, and
   report calibrated FDR (roadmap item 2)
4. Commit the results so regressions are visible in review
