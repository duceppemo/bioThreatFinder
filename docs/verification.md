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

Current state: **286 passed, 0 failed, 0 skipped** (2026-10-20) with `samtools` on
PATH and `BTF_TAXDUMP` set. The bare `python3 tests/test_pipeline.py` reports
**267 passed, 3 skipped** - the `coverage` and `statecols` checks need `samtools`, and
the `screen`/`leadjoin` writers need a taxdump. The per-script `--self-test` suites run
as single checks from here, so `selftest:build_marker_db` alone covers 60+ assertions,
including D30's eleven on the protein query ladder.

## What the tests assert

| Group | Checks | The behaviour being pinned down |
|---|---|---|
| `selftest` | 5 | The self-tests inside `btf_yaml`, `marker_rules`, `evaluate_benchmark`, `build_marker_db`, `mask_pack_refs` |
| `tier` | 9 | All six decision paths: true positive, conserved-repeat artifact → `attributed_to_neighbour`, below-floor lead, control-suppressed, control's own background, and a database gap surfacing as `no_reference` rather than a negative |
| `marker` | 10 | **Identical** alignment evidence yielding four different tiers on marker state alone: both plasmids → confirmed; pXO2 confidently absent → held; pXO2 **untested** → held *and named as a database gap*; *B. cereus* cereulide present → veto. Plus BCG → `vaccine_strain`, and disclosure when a veto could not be evaluated |
| `d17a` | 7 | A taxon with no marker rule block is capped at `supported`, and the escape hatch the cap prints actually works (D37 - it did not: `RULE_KEYS` never contained `no_markers_available`, and the test named for it wrote a rules file it never passed to the scorer). A declared gap lifts the cap, still discloses that the marker tier never ran, is NOT lifted by `false`, and refuses a mistyped value rather than reading it as permission |
| `wright` | 6 | The Wright et al. composite: genuine call confirms; **96.2% identity held** (genus right, species wrong); 18% mapped fraction held; and a genuinely divergent CCHFV at 89.2% still confirming via the virus-specific floor |
| `coverage` | 19 | Nine alignment metrics exact against a hand-built BAM with known flags, MAPQs, strands and duplicates — plus masking shrinking the breadth denominator and excluding masked-region reads while still reporting them |
| `screen` | 17 | Multi-classifier parsing; a piled-up taxon flagged not dropped; HLL cardinality read from the clade rather than summed; unclassified fraction reported and flagged. Plus the four stage-1 coverage states (D33): an arm nobody asked for is `not_run` and never a failure, an **enabled** arm that produced nothing is `no_output`, an arm with no read ledger reports blank rather than `0`, a run where every arm died still writes its table, and a misspelled `--expected-arms` value is rejected rather than becoming an unenforced expectation. Plus D36: an arm that searched only the PACK answers a clean sample with an empty report, and that is the negative answer rather than a failure - while the same arm's *absent* file is still `ARM_PRODUCED_NOTHING` |
| `pack` | 18 | Rules compile and reference only real markers; scenario file parses; all seven NYSDOH categories present; plumbing positive control exists; every scenario taxon has a genome source; identity floors are group-aware; **every retrievable marker can build at least one query** (D30 — `henipa_cedar` built none, and a row that is never searched must not be filed as a retrieval failure) |
| `aacov` | 8 | Protein coverage is per allele, not over their summed length (D31). A full-length hit on one of four alleles is PRESENT at coverage 1.0; the same hit at 20%, and 30% on each of four alleles, are both ABSENT; a target missing from the FASTA is warned about rather than silently scoring 0 |
| `depth` | 13 | Marker absence is only called when depth could have shown it; below that it is UNKNOWN, not ABSENT |
| `eval` | 9 | The benchmark evaluator's tier→outcome mapping, LoD derivation and expectation checking. Its self-test also now pins that a sample the run never analysed is EXCLUDED from grading rather than scored as a miss, that one which ran and produced nothing is still a miss, and that the old behaviour is kept when no samplesheet is given (D35 — graded against the whole manifest, this run reported a false FAIL on the one notifiable organism it had called correctly). D38 added ten more: a `max_tier` specificity cap is breached only by a tier that ASSERTS PRESENCE, so a neighbour attribution, a marker veto and a group lead each PASS one and each say why, while `supported`, `confirmed_in_silico`, `indistinguishable_complex` and `vaccine_strain` still FAIL — the rank-only comparison failed the pipeline's strongest exonerations, `marker_veto` outranking a POSITIVE tier |
| `absent` | 6 | `truth_absent` is the full watchlist minus what is genuinely present, so a false positive on an unlisted organism cannot be scored `unlisted` and dropped from specificity (D14) |
| `platform` | 3 | Every finding carries ITS OWN row's platform; a short-circuiting row cannot inherit the previous one's |
| `shared_assembly` | 7 | Taxa sharing one NCBI assembly carry identical masking numbers across all eight columns, and neither is the vacuous `genome_bp=0` row. Each partner now also NAMES the other in `SHARED_SEQUENCE` (D34 - the flag is derived from the manifest's `accession`, so it no longer depends on how the script was invoked). Plus the two invariants whose violation was invisible: **every pack reference has a masking result**, and no reference is lost to taxid-keyed bookkeeping - this test itself read the manifest into `{taxid: row}` and so discarded eleven of taxid 11266's twelve rows |
| `workflow` | 6 | Every script the workflow invokes is in `bin/` and executable; no shared placeholder for optional inputs; stage 2 does not re-filter on assigned reads |
| `report` | 13 | **The first coverage `make_report.py` has ever had** - which is part of why an entire input could go missing from it unnoticed (D33). A stage-1 arm that was enabled and produced nothing is bannered above the findings, not buried; the unclassified fraction reaches the page; bracken's by-design exit carries no alarm; and a report rendered without a screen summary omits the section rather than drawing an empty one that reads as all-clear. D36 added the run-level view that pays for that suppression: an arm empty on EVERY sample it searched is bannered as a database problem, because an empty or truncated database loads without error and is silent per-sample |
| `ssrole` | 12 | The samplesheet declares which samples are controls, so it DRIVES the subtraction and not merely the count (D39). A sheet-declared control suppresses even when the pipeline's own rows call it a sample - the re-score shape, where declaring a control used to assert *"all clean"* on every finding while reporting that control as a `supported` detection. The sheet is authority in both directions, a declared control that never ran is not credited as a clean one, and with no read ledger to tell a spotless NTC from an absent one the question is left unasked |
| `markergene` | 12 | A marker's sequences must be THE GENE THE CATALOGUE NAMES, and three shipped markers were not (D40). The nucleotide and protein tiers are retrieved by separate ladders, so they can land on different genes while the manifest says `built` and marker_evidence treats the marker as fully testable - `min_coverage` then being a fraction of the wrong reference. The matcher no longer lets a one-letter CDS name claim a word alias (`[gene=P]` claimed `nipah_L` because *"polymerase".startswith("p")*, so the L-polymerase marker's nt reference is 2130 bp of the Nipah P gene, in a BSL-4 agent's `confirm_requires`), an exact gene match now outranks an operon-partner one (`ypest_pmt1_caf1` holds caf1A/M/R and not caf1), and the permissiveness that is deliberate - `capBCADE` for `capB`, `ORF1ab` for `ORF1a`, `pag` for `pagA` - still works. Plus the invariant across all 115 checkable references, with a SELF-RETIRING pending list: the test also asserts each pending defect is still real, so a rebuild forces the entry out |
| `salm` | 4 | The Salmonella block is species confirmation, not pathotype (D40). All three loci present confirms and the finding names the rule; dropping `ssaV` does NOT confirm, because invA and ttrC are both genus-level and S. bongori carries both; and an UNTESTED `ssaV` holds it just as firmly as an absent one. The fixture carries a read ledger on purpose - without it the alignment-length floor falls back to a flat 100 bp and a 2x101 library's 99 bp mean fails the very gate D37 found, so the test would have passed for the wrong reason |
| `statecols` | 3 | The states file's columns must not depend on which marker sorts first (D40). `marker_evidence.py` took its header from `results[0]`, and a record only carries `nt_depth`, `nt_expected_breadth`, `nt_n_seqs` and `nt_best_ref` if it reached the nucleotide measurement block - so a catalogue whose first row is UNKNOWN produced a header without them and `extrasaction="ignore"` dropped them from every row, silently, including rows that had measured them. `nt_expected_breadth` is the input to the depth gate, the rule that stops a false ABSENT |

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
  As of 2026-10-16 all **371** pack references have a result: the twelve Filoviridae
  group-entry files had never been masked (D34) and came back at 0 bp. The script now
  reconciles the whole pack and exits non-zero when a reference has no result, so the
  count no longer has to be checked by hand.
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
  `scripts/fetch_mock_communities.py` exist for that.
- **(superseded) Real-data accuracy IS now measured, on ZymoBIOMICS D6300.**
  `runs/zymo_rl100_rescore` (2026-10-20, 11 Illumina libraries): **no false
  confirmations on absent taxa**, and tier PPV is monotone - `supported` 0.600
  (5 calls), `attributed_to_neighbour` 0.111, `screened_only` 0.000.
  **`confirmed_in_silico` has zero calls on real data**, so the question the PPV table
  exists to answer is still unanswered: the top tier is unreachable for the 32 of 90
  watchlist taxa with no marker rule block (D37), and every other `supported` call is
  held by the 0.98 illumina identity floor, which a real *L. monocytogenes* isolate
  misses at 0.9643 through strain divergence from the pack reference (D38). Both
  figures are a FLOOR, not a measurement, because **no real-data samplesheet declares
  a negative control** - so two `supported` calls that are almost certainly library
  cross-contamination cannot be told from false positives. See D37, D38.
- **`confirmed_in_silico` HAS now fired on real data, on curated marker evidence**
  (D40, 2026-10-20). *S. enterica* in the `ERR2935848` isolate library passes all six
  alignment gates and all three markers of a real `marker_rules.yaml` block -
  `salm_invA` + `salm_ttrC` + `salm_ssaV` at breadth 1.000, identity 0.985-0.989 and
  475-605x depth. Exactly one row moved across all 82 findings. D37 reached this tier
  once before on a SCRATCH `no_markers_available` declaration that was never committed;
  this is the first time on curated evidence. The call still discloses the missing
  negative control and the missing KrakenUniq arm. NOT a measurement of accuracy: it is
  one organism in its own isolate library, and 31 of 90 taxa still have no block.
- **What ONE declared control is worth has now been measured** (D39, 2026-10-20), on
  the same 11 libraries and the same alignment evidence: `supported` PPV **0.600 ->
  0.750**, false positives **2 -> 1**, false confirmations 0, and **every per-scenario
  detection rate unchanged** - all 3 suppressed rows were truth-absent, so the control
  cost no sensitivity on this data. The remaining false positive is Listeria in the
  B. subtilis library, surviving on a **5.11x** margin against the 5.0x requirement
  while carrying Listeria's own `hly` and `prfA` markers: real contaminating DNA, not
  spurious mapping. This is a MEASUREMENT OF THE INPUTS, not of a committed change -
  the yeast library is another organism's library from the same run, not a no-template
  control, and no samplesheet in the repo was edited. Getting here first required
  fixing the scorer: declaring a control on the re-score path silenced the
  missing-control warning and subtracted nothing, putting *"1 negative control(s) in
  run, all clean"* on all 82 findings of a run whose declared control was itself
  reported as a `supported` L. monocytogenes detection at 355.5 RPM.
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
