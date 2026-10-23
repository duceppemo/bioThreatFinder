# Verification

What has actually been checked, how to re-run it, and — the part that matters most —
**what has not been checked at all.**

## Run it

```bash
python3 tests/test_pipeline.py          # ~30 s, no databases needed
python3 tests/test_pipeline.py --list
python3 tests/test_pipeline.py -k marker
```

Needs `samtools` for the coverage cases (skipped with a message if absent) and `$BTF_DB`
pointing at built databases for the cases that compare a shipped artifact against the repo
— the taxdump, the refs manifest, the sketch, the marker database and the Kraken2
taxonomy. All of those skip with a message naming what to set. Everything else is
fixture-driven.

Current state: **477 passed, 0 failed, 0 skipped** (2026-10-23) with `samtools` on PATH
and `$BTF_DB` set. On a bare clone with neither, `python3 tests/test_pipeline.py` reports
**416 passed, 0 failed, 14 skipped** - 3 of those skips want `samtools`, the other 11 want
built databases. Per-case overrides (`BTF_TAXDUMP`, `BTF_REFS_DIR`, `BTF_MARKER_DB`,
`BTF_REFS_MANIFEST`) still point one case at a copy elsewhere. The per-script `--self-test` suites run
as single checks from here, so `selftest:build_marker_db` alone covers 60+ assertions,
including D30's eleven on the protein query ladder.

## What the tests assert

| Group | Checks | The behaviour being pinned down |
|---|---|---|
| `selftest` | 6 | The self-tests inside `btf_yaml`, `marker_rules`, `evaluate_benchmark`, `build_marker_db`, `mask_pack_refs`, `measure_reference_span`. The last covers the four ways D43's measurement could print a confident table instead of an error: the plain shipped filename may be claimed only by a taxid's FIRST reference, the cover is seeded and prefers single-sequence assemblies on ties, a near-zero shared-hash count is counted as unmeasurable rather than reported as a precise low ANI, and every external binary is RUN before use (D41) |
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
| `mref` | 8 | Breadth dilution follows from having SEVERAL references, not from rank (D42). A species with two single-sequence assemblies is judged on the one its reads are actually on, and the substitution changes the TIER rather than only printing a note. A single-reference taxon gets no note at all. A DRAFT assembly blocks the substitution instead of letting one contig stand for a genome, and says the row may be under-called. With no refs manifest the question is left unasked |
| `stale` | 13 | Existence is not freshness (D45). A competitive bundle is rebuilt when its INPUTS changed, not merely when it is absent: a sidecar records the assemblies actually concatenated, keyed on accession, so a reference added, removed or swapped is caught with every mtime settled - the shape D42 hit, which mtime alone reports as current - while mtime still catches content changed under a known name. Missing and EMPTY are both stale. Indexes likewise: an index older than its bundle is the same disagreement one step on. Both decisive cases were run against the old existence-only predicate, which calls each of them current |
| `sketch` | 19 | D45's predicate arriving at the one artifact that had no check: the sylph sketch (D46). The shell `--check` is run against a STUB `sylph inspect`, so every verdict is exercised offline - a reference never sketched is caught with all mtimes settled (the state that actually shipped for 22 days), one newer than the sketch is caught by mtime, a sketched genome no longer in `refs/` is caught, and a database at the wrong `c` is stale too because it matches nothing while loading cleanly. Undetermined stays distinguishable from a verdict: an unreadable database exits 1, an absent one exits 3. The guard is asserted to be CALLED in the sylph branch, not merely defined; and no guard may count a list it deduplicated with Groovy's in-place `.unique()`, which made both this guard and the masking guard report "1 reference" for the twelve-file state they exist to report. A build RECORDS the genome set beside the database and `--stamp` records what the DATABASE holds rather than what `refs/` holds - derived from `refs/` it would agree by construction, the one thing it must be able to disagree with - and a record disagreeing with its database is caught, because a run trusts that file and cannot look inside the `.syldb`. The builder and the run are asserted to name the SAME record file: if those drift nothing fails, the builder writes a record nobody reads, and every run falls back to the timestamps-only check the record replaced. Finally the shipped artifacts: no reference newer than the sketch, its recorded genome count equal to the pack's reference count, and the record naming EXACTLY the pack's 375 references |
| `markerdb` | 23 | The same predicate at the marker database, whose input is `markers.tsv` (D47), and the mtime/set relationship is INVERTED from the sketch's: here the sets are blind to a marker **edited** under an id the manifest already knows - every set agrees while the sequence underneath is the one the last build retrieved, which is D41's wrong-gene shape - so mtime is the only signal that sees it, and the case is run with the mtime comparison removed, where it reports FRESH. A claimed tier the combined artifact does not hold is caught because that one fails UNSAFE: `marker_evidence.py` believes the manifest, searches, finds nothing and reports ABSENT, asserting the marker was looked for. Presence is compared and never sequence COUNTS - `concat_files()` de-duplicates within a marker while `n_nt_seqs` counts the per-marker archive, and the count-comparing version was run: it calls the correct shipped database stale on **116 of 166** markers, which is how a guard gets switched off. An index older than its source is caught, a missing one too; an absent database is stale (3) while an absent catalogue is undetermined (1), so a verdict stays distinguishable from a failure. `--check` is asserted to answer BEFORE `load_rankedlineage()`, because a read-only question costing minutes and gigabytes is one nobody asks twice. The guard is asserted to be CALLED in the marker branch, and the column parser to be a **top-level `def`**: written first as a closure in a local `def` and called as `idsOf(...)`, it failed to COMPILE under Nextflow 26.04's strict syntax - every run with markers enabled, which is the default - and that was found by extracting the function and running it, not by reading it. Finally the shipped database itself, checked in full on every commit, which is where the two build-integrity comparisons the launch guard leaves to `--check` are actually exercised |
| `envs` | 7 | The three build environments as FILES rather than a paragraph (they were documented only in prose, so an outside user had to reconstruct them from a description). These are curated specs, not `conda env export` output, and the invariant is why: the environment they were derived from still holds **kraken2 2.0.8-beta**, because installing kraken2 beside mmseqs2/hmmer/mafft makes the solver pick it - and it has no `--report-minimizer-data`, so an index built with it works while carrying none of the evidence stage 1 scores on. An export would therefore reproduce the trap. `btf-refs` is asserted NOT to declare kraken2 or bracken; both pins are asserted load-bearing rather than hygiene (kraken2 must match the version that BUILT the index or the format may not read back, and sylph's pin must equal the one in `modules/local/detect.nf` because a sketch must be readable by the version that queries it); and every file is asserted to declare only conda-forge and bioconda, since `defaults` imposes licence conditions on a reader and the originals had inherited unrelated channels from one developer's global config |
| `idxcover` | 6 | The Kraken2 index's question is COVERAGE, not freshness (D48): a watchlist taxon with no node in its pruned taxonomy can never appear in its report, and the arm still records `ok`. **19 of 90 are absent** on the shipped index - helminths, several protozoa, most dimorphic fungi - which follows from PlusPF's library selection and not from staleness, so a rebuild changes nothing. The number is pinned in BOTH directions: more is not automatically a regression and fewer is not automatically a fix, so either way the instruction is to re-measure and update D48 rather than adjust the number until the test passes. The bound is guarded harder than the gap: all 19 are sketched for sylph, which is on by default with `combining_rule: any_classifier_or`, so one arm is blind rather than the stage - and a taxon with NO node and NO sketched reference would be invisible to stage 1 entirely, which the check reports separately and this test fails on. Two overclaims are blocked: the pruned taxonomy (62,771 nodes) is never compared to the pinned taxdump's 3,016,439 - that comparison calls every correct index broken - and a present node is stated as necessary, not sufficient, since it can exist because a relative carried the genome. `taxo.k2d` is read with no kraken2 on the head node, and the layout arithmetic must add up to the byte or the reader REFUSES: a wrong stride does not fail, it reads the wrong 8 bytes of every record and reports them as taxids confidently |
| `extrarefs` | 11 | Extra references come from a MEASURED cover, not from curation rank (D44). Both input formats parse; a missing cover file is FATAL rather than silently empty, because loading nothing would build a one-reference-per-taxon pack while the operator believed a cover had been applied; the curated assembly stays first and keeps the plain filename earlier measurements used; pins follow in the cover's order and are not re-sorted, which would discard the measurement; and a pin already in the ranked list is not held twice - the two-rows-one-accession failure of D42, reachable again through this path, and verified to fail against the naive implementation |
| `refs` | 3 | Invariants over the shipped reference pack, all three violated on 2026-10-21 and each verified to fail against a reproduction of the broken state (D42): a `(taxid, accession)` pair appears ONCE; a taxon's row count equals its distinct accession count - `--max-per-taxon N` deduped candidates by FILENAME, so it re-downloaded accessions already held and credited a taxon with N references while holding one genome N times; and no competitive bundle is older than the genomes it is built from, because the FASTA and indexes are skipped whenever the output merely exists, leaving the manifest saying two references and the aligner using one |

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

- **Watchlist resolution** against a live taxdump: 90/90 (89/89 when first run, before the
  `Filoviridae` group entry), while flagging 18 changed virus ranks, 11 ICTV renames
  resolving via equivalent names, 4 malformed neighbour entries.
- **Competitive bundles**: all 90 scanned for identical sequence under two taxid prefixes
  — 1 found (*F. tularensis*/*F. novicida*), now 0. The check is built into
  `build_pack_refs.py` as a post-write step.
- **Taxonomy drift**: the previously-default snapshot measured 201 days old, +288,764
  taxids behind current, 2,581 descendant nodes short inside watchlist clades.
- **Masking**: 20,040,354 bp over 375 references (0.348%); 104.2 kb vector across 45
  contigs. As of 2026-10-21 all **375** pack references that have sequence carry a result;
  the twelve Filoviridae group-entry files had never been masked (D34, when the pack held
  371) and came back at 0 bp. The script now
  reconciles the whole pack and exits non-zero when a reference has no result, so the
  count no longer has to be checked by hand.
- **The sylph pack sketch was reconciled against the refs manifest** (2026-10-21) and was
  **16 references short** — all twelve `11266_Filoviridae*` files and four *Listeria*,
  added 2026-10-14 and 2026-10-21 against a sketch built 2026-09-29. sylph is on by
  default, so for 22 days that arm could not raise the Filoviridae group entry whose whole
  purpose is catching a divergent filovirus, and nothing said so: the arm ran, produced a
  usable report, and `screen_summary.tsv` recorded it `ok`. Rebuilt and verified by
  `sylph inspect` on both sketches and diffing the genome sets — 359 → 375, exactly the 16
  expected files added, none removed, c/k/min-spacing 200/31/30 matching the query side.
  ~~**Not automated:**~~ **Now automated (D46)**, and the automation was verified against
  this same artifact: `scripts/build_sylph_db.sh --check` pointed at the retained
  `human_pack_refs.bak-2026-09-29.syldb` names the 16 missing references and exits 3, while
  the live sketch reports FRESH. The same comparison now happens **at launch**: the sketch
  carries a record of its genome set (`<db>.genomes.txt`, read back from `sylph inspect`)
  and `warnIfSketchStale()` checks it against the pack with no sylph needed. Verified by
  extracting the shipped function and running it under Nextflow against the real manifest
  with the old sketch's own record — it names all 16 at launch — and against a synthetic
  pack where the database is newer than every reference, so mtime reports current and only
  the record catches the unsketched one.
- **The reference pack's reach across each target species was measured without a
  specimen** (2026-10-21, D43). 13,582 known assemblies of the 26 bacterial targets,
  every one downloaded, each measured against the masked reference the pipeline actually
  aligns to. 17 of 26 species already span their own diversity on one reference; 9 do
  not. The verdict comes from simulating reads through the pipeline's own aligner,
  bundle and `coverage_stats.py`, NOT from an ANI cut-off, because the ANI at which
  `mean_identity` crosses 0.98 is species-dependent (~0.973 to ~0.981 over the three
  species calibrated). Instrument controls: each shipped reference scores ANI 1.0000
  against itself, and a genome simulated against itself returns `mean_identity` 0.9980,
  so the ceiling sits under two points above the floor. Report: `docs/reference_span.tsv`.

- **A second reference for one species was measured, not assumed** (2026-10-21, D42).
  Same reads against a one- and a two-reference bundle: `mean_identity` **0.9641 ->
  0.9947**, `best_reference_breadth` **0.877 -> 0.982**, and **6.2 M more reads map at
  all, a 48% gain** - they could not previously align to the wrong-lineage genome. The
  0.9641 control reproduces the production run's 0.9643, which is what makes the
  comparison valid. The 0.98 species-discrimination floor was NOT lowered and no
  per-taxon override was added. Raw breadth falls 0.877 -> 0.544 because the denominator
  now spans both genomes, which is what the `mref` substitution exists to handle.
- **The marker database was rebuilt end to end** (2026-10-21, D41): 166 rows forced,
  143 -> 146 `built`, zero status regressions, 33 profiles unchanged. Verified not from
  the manifest alone but by realigning a real library against the rebuilt references,
  deriving states from the production manifest and re-scoring - `confirmed_in_silico`
  fires on the database the pipeline actually loads. The first attempt was ROLLED BACK:
  a local `hmmbuild` that resolved but could not run degraded 33 markers to `partial`.
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
  sit under `$BTF_DB/runs/` or `$BTF_DB/simulated/*/bakeoff/`. Those runs are
  where D12, D13, D14 and D15 came from, and where the Toxoplasma false positive was found
  and then fixed: `lod_v4` is the first batch carrying a host-matrix control, and its
  `supported` tier is empty where `lod_v3` had 18 *T. gondii* calls.
- **No threshold is measured.** Every number in `packs/*/pack.yaml` is a literature
  default (Wright et al., UKHSA mSCAPE, Oskolkov) or a reasoned guess.
- **Conda environments HAVE been exercised.** 13 `.nextflow-conda/env-*` directories
  exist and every run above used `-profile mamba`. `docs/operations.md` records all 13 specs
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
