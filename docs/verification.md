# Verification

What has actually been checked, how to re-run it, and — the part that matters most —
**what has not been checked at all.**

## Run it

```bash
python3 tests/test_pipeline.py          # ~10 s, no databases needed
python3 tests/test_pipeline.py --list
python3 tests/test_pipeline.py -k marker
```

Needs `samtools` for the coverage cases (skipped with a message if absent) and a taxdump
for the screen cases (`BTF_TAXDUMP=/path/to/taxdump`, skipped otherwise). Everything else
is fixture-driven.

Current state: **76 checks, all passing.**

## What the tests assert

| Group | Checks | The behaviour being pinned down |
|---|---|---|
| `selftests` | 4 | The self-tests inside `btf_yaml`, `marker_rules`, `evaluate_benchmark`, `build_marker_db` |
| `tiering` | 6 | All six decision paths: true positive, conserved-repeat artifact → `attributed_to_neighbour`, below-floor lead, control-suppressed, control's own background, and a database gap surfacing as `no_reference` rather than a negative |
| `marker` | 7 | **Identical** alignment evidence yielding four different tiers on marker state alone: both plasmids → confirmed; pXO2 confidently absent → held; pXO2 **untested** → held *and named as a database gap*; *B. cereus* cereulide present → veto. Plus BCG → `vaccine_strain`, and disclosure when a veto could not be evaluated |
| `wright` | 6 | The Wright et al. composite: genuine call confirms; **96.2% identity held** (genus right, species wrong); 18% mapped fraction held; and a genuinely divergent CCHFV at 89.2% still confirming via the virus-specific floor |
| `coverage` | 16 | Nine alignment metrics exact against a hand-built BAM with known flags, MAPQs, strands and duplicates — plus masking shrinking the breadth denominator and excluding masked-region reads while still reporting them |
| `screen` | 4 | Multi-classifier parsing; a piled-up taxon flagged not dropped; HLL cardinality read from the clade rather than summed; unclassified fraction reported and flagged |
| `pack` | 7 | Rules compile and reference only real markers; scenario file parses; all seven NYSDOH categories present; plumbing positive control exists; every scenario taxon has a genome source; identity floors are group-aware |

## It has already caught a real bug

On its first run, `marker:BCG` failed. Adding the veto/vaccine CAVEAT disclosures had
broken an `elif` chain, so the `confirm_requires` branch **overwrote** a `vaccine_strain`
determination:

```
- VACCINE/ATTENUATED STRAIN indicated ...        <- tier set correctly
- HELD at supported: required markers not found  <- then overwritten
```

A BCG vaccine strain would have been reported as `supported` — a possible tuberculosis
detection — instead of a stand-down. Wrong direction for a clinical call, invisible in
normal use, and found because the behaviour was written down as an executable assertion
instead of a commit message.

Two of the three initial failures were bugs in the *tests*, both worth recording: findings
were keyed by sample alone, but a `no_reference` row shares the sample name; and the
artifact fixture sat below the read floor after the UKHSA-derived recalibration, so it
never reached the gate it was meant to exercise.

## Verified by hand, outside the suite

One-offs against real data, not automated:

- **Watchlist resolution** against a live taxdump: 89/89, while flagging 18 changed virus
  ranks, 11 ICTV renames resolving via equivalent names, 4 malformed neighbour entries.
- **Competitive bundles**: all 89 scanned for identical sequence under two taxid prefixes
  — 1 found (*F. tularensis*/*F. novicida*), now 0. The check is built into
  `build_pack_refs.py` as a post-write step.
- **Taxonomy drift**: the previously-default snapshot measured 201 days old, +288,764
  taxids behind current, 2,581 descendant nodes short inside watchlist clades.
- **Masking**: 20,028,976 bp over 359 genomes; 93.9 kb vector across 44 contigs.
- **Guards refuse correctly**: GTDB-namespaced Metabuli index rejected; wrong refs base
  rejected; 201-day taxdump rejected by both the resolver and the Kraken2 builder;
  kraken2 2.0.8-beta rejected for lacking `--report-minimizer-data`.

## What is NOT verified — read this part

- **No accuracy measurement of any kind.** No ground-truth data has been through this
  pipeline. Nothing here says the design is *right*, only that it behaves as designed.
- **The pipeline has never run end to end.** Every process is individually plausible and
  the DAG parses; no sample has been processed.
- **No threshold is measured.** Every number in `packs/*/pack.yaml` is a literature
  default (Wright et al., UKHSA mSCAPE, Oskolkov) or a reasoned guess.
- **No conda environment resolved by the pipeline has been exercised.** The pins in
  `modules/local/*.nf` were corrected but never installed; several are major-version jumps.
- **Sensitivity is untested in both directions.** No limit of detection, and no false-negative
  rate — the tests confirm that things are *held back* correctly, which is half the story.
- **Simulated data would flatter us anyway.** See `docs/benchmarking.md` for what
  simulation cannot tell you.

## The next real test

`scripts/simulate_benchmark.py` → run the pipeline → `bin/evaluate_benchmark.py`. That
produces the first numbers about accuracy rather than behaviour. Until it runs, this
document is the honest ceiling on what can be claimed.
