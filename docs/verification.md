# Verification

What has actually been checked, how to re-run it, and — the part that matters most —
**what has not been checked at all.**

## Run it

```bash
python3 tests/test_pipeline.py          # ~30 s, no databases needed
python3 tests/test_pipeline.py --list
python3 tests/test_pipeline.py -k marker
```

Needs `samtools` for the coverage cases (skipped with a message if absent) and a taxdump
for the screen cases (`BTF_TAXDUMP=/path/to/taxdump`, skipped otherwise). Everything else
is fixture-driven.

Current state: **203 passed, 0 failed, 2 skipped** (2026-10-15). The 2 skips are the
`coverage` checks that need `samtools` on PATH. The per-script `--self-test` suites run
as single checks from here, so `selftest:build_marker_db` alone covers 60+ assertions,
including D30's eleven on the protein query ladder.

## What the tests assert

| Group | Checks | The behaviour being pinned down |
|---|---|---|
| `selftest` | 4 | The self-tests inside `btf_yaml`, `marker_rules`, `evaluate_benchmark`, `build_marker_db` |
| `tier` | 9 | All six decision paths: true positive, conserved-repeat artifact → `attributed_to_neighbour`, below-floor lead, control-suppressed, control's own background, and a database gap surfacing as `no_reference` rather than a negative |
| `marker` | 10 | **Identical** alignment evidence yielding four different tiers on marker state alone: both plasmids → confirmed; pXO2 confidently absent → held; pXO2 **untested** → held *and named as a database gap*; *B. cereus* cereulide present → veto. Plus BCG → `vaccine_strain`, and disclosure when a veto could not be evaluated |
| `wright` | 6 | The Wright et al. composite: genuine call confirms; **96.2% identity held** (genus right, species wrong); 18% mapped fraction held; and a genuinely divergent CCHFV at 89.2% still confirming via the virus-specific floor |
| `coverage` | 19 | Nine alignment metrics exact against a hand-built BAM with known flags, MAPQs, strands and duplicates — plus masking shrinking the breadth denominator and excluding masked-region reads while still reporting them |
| `screen` | 12 | Multi-classifier parsing; a piled-up taxon flagged not dropped; HLL cardinality read from the clade rather than summed; unclassified fraction reported and flagged. Plus the four stage-1 coverage states (D33): an arm nobody asked for is `not_run` and never a failure, an **enabled** arm that produced nothing is `no_output`, an arm with no read ledger reports blank rather than `0`, a run where every arm died still writes its table, and a misspelled `--expected-arms` value is rejected rather than becoming an unenforced expectation |
| `pack` | 18 | Rules compile and reference only real markers; scenario file parses; all seven NYSDOH categories present; plumbing positive control exists; every scenario taxon has a genome source; identity floors are group-aware; **every retrievable marker can build at least one query** (D30 — `henipa_cedar` built none, and a row that is never searched must not be filed as a retrieval failure) |
| `aacov` | 8 | Protein coverage is per allele, not over their summed length (D31). A full-length hit on one of four alleles is PRESENT at coverage 1.0; the same hit at 20%, and 30% on each of four alleles, are both ABSENT; a target missing from the FASTA is warned about rather than silently scoring 0 |
| `depth` | 13 | Marker absence is only called when depth could have shown it; below that it is UNKNOWN, not ABSENT |
| `eval` | 9 | The benchmark evaluator's tier→outcome mapping, LoD derivation and expectation checking |
| `absent` | 6 | `truth_absent` is the full watchlist minus what is genuinely present, so a false positive on an unlisted organism cannot be scored `unlisted` and dropped from specificity (D14) |
| `platform` | 3 | Every finding carries ITS OWN row's platform; a short-circuiting row cannot inherit the previous one's |
| `shared_assembly` | 4 | Taxa sharing one NCBI assembly carry identical masking numbers across all eight columns, and neither is the vacuous `genome_bp=0` row |
| `workflow` | 6 | Every script the workflow invokes is in `bin/` and executable; no shared placeholder for optional inputs; stage 2 does not re-filter on assigned reads |
| `report` | 7 | **The first coverage `make_report.py` has ever had** - which is part of why an entire input could go missing from it unnoticed (D33). A stage-1 arm that was enabled and produced nothing is bannered above the findings, not buried; the unclassified fraction reaches the page; bracken's by-design exit carries no alarm; and a report rendered without a screen summary omits the section rather than drawing an empty one that reads as all-clear |

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
- **Masking**: 20,040,144 bp over 359 genomes (0.349%); 104.4 kb vector across 45 contigs.
- **Guards refuse correctly**: GTDB-namespaced Metabuli index rejected; wrong refs base
  rejected; 201-day taxdump rejected by both the resolver and the Kraken2 builder;
  kraken2 2.0.8-beta rejected for lacking `--report-minimizer-data`.

## What is NOT verified — read this part

- **Accuracy IS now measured, on simulated data only.** Limit of detection for
  *Y. pestis* on Illumina is **1000 RPM**, a cliff rather than a gradient (0.000 detection
  at 1/10/100, 1.000 at 1000/10000), measured over `runs/lod_v1` and `runs/lod_v2_kraken2`
  against the `lod_v1` manifest. sylph alone and Kraken2+sylph give the SAME limit. See
  `docs/validation.md` and decisions D12/D13/D14. What remains unmeasured is accuracy on
  REAL sequencing data - `benchmarks/scenarios/real_mock_communities.yaml` and
  `scripts/fetch_mock_communities.py` exist for that and nothing has been downloaded.
- **The pipeline HAS run end to end, repeatedly** - `plumbing_001`, `plumbing_v2`,
  `lod_v1`, `lod_v2_kraken2`, `mmap_test`, `hostqc_test`, `metabuli_smoke`,
  `lod_v3_threearm`, `lod_v4`, the `divergent_v1` bake-off and (in progress) `d12fine` all
  sit under `/media/8tb/bioThreatFinder/runs/` or `simulated/*/bakeoff/`. Those runs are
  where D12, D13, D14 and D15 came from, and where the Toxoplasma false positive was found
  and then fixed: `lod_v4` is the first batch carrying a host-matrix control, and its
  `supported` tier is empty where `lod_v3` had 18 *T. gondii* calls.
- **No threshold is measured.** Every number in `packs/*/pack.yaml` is a literature
  default (Wright et al., UKHSA mSCAPE, Oskolkov) or a reasoned guess.
- **Conda environments HAVE been exercised.** 13 `.nextflow-conda/env-*` directories
  exist and every run above used `-profile mamba`. `docs/STATUS.md` records all 13 specs
  put through the solver: 12 solved, samtools 1.24 did not, and samtools is now pinned
  1.22. The original note follows for the pins it still describes.
- **(historical) No conda environment resolved by the pipeline has been exercised.** The pins in
  `modules/local/*.nf` were corrected but never installed; several are major-version jumps.
- **(superseded - an LoD now exists, see above)** Sensitivity is untested in both directions. No limit of detection, and no false-negative
  rate — the tests confirm that things are *held back* correctly, which is half the story.
- **Simulated data would flatter us anyway.** See `docs/benchmarking.md` for what
  simulation cannot tell you.

## The next real test

`scripts/simulate_benchmark.py` → run the pipeline → `bin/evaluate_benchmark.py`. That
produces the first numbers about accuracy rather than behaviour. Until it runs, this
document is the honest ceiling on what can be claimed.
