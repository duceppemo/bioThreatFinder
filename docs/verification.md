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

Current state: **517 passed, 0 failed, 0 skipped** (2026-10-26) with `samtools` on PATH
and `$BTF_DB` set. On a bare clone with neither, `python3 tests/test_pipeline.py` reports
**451 passed, 0 failed, 16 skipped** - 3 of those skips want `samtools`, the other 13 want
built databases. Per-case overrides (`BTF_TAXDUMP`, `BTF_REFS_DIR`, `BTF_MARKER_DB`,
`BTF_REFS_MANIFEST`) still point one case at a copy elsewhere. The per-script `--self-test` suites run
as single checks from here, so `selftest:build_marker_db` alone covers 60+ assertions,
including D30's eleven on the protein query ladder.

## What the tests assert

| Group | Checks | The behaviour being pinned down |
|---|---|---|
| `selftests` | 9 | The self-tests inside `btf_yaml`, `marker_rules`, `evaluate_benchmark`, `build_marker_db`, `mask_pack_refs`, `measure_reference_span`, `check_index_coverage`, `host_depletion_stats` and `to_pavian`. The last two were added by D62's review: `host_depletion_stats` had a passing `--self-test` that nothing ever invoked, and `to_pavian.py` was the only pipeline-invoked script with no coverage at all — which is the wrong script to leave uncovered, because it is written to DEGRADE rather than die, so a converter handed the wrong file produced fabricated taxids under exit 0. `measure_reference_span` covers the four ways D43's measurement could print a confident table instead of an error: the plain shipped filename may be claimed only by a taxid's FIRST reference, the cover is seeded and prefers single-sequence assemblies on ties, a near-zero shared-hash count is counted as unmeasurable rather than reported as a precise low ANI, and every external binary is RUN before use (D41) |
| `tiering` | 9 | All six decision paths: true positive, conserved-repeat artifact → `attributed_to_neighbour`, below-floor lead, control-suppressed, control's own background, and a database gap surfacing as `no_reference` rather than a negative |
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
| `findcols` | 3 | The same defect as `statecols`, in `score_evidence.py` (D62). The findings header came from `findings[0]`, and the `no_reference` rows - a watchlist taxon whose reference was never built - are assembled from a different literal that omitted `best_neighbour_uniformity`, so a run whose findings are ALL database gaps wrote 36 columns instead of 37 and dropped the neighbour-discrimination measure D18 put in the TSV so the gate could be audited from it. A mixed run was safe only by accident: `no_reference` is the lowest tier and sorts last, so `findings[0]` happened to be complete. **The union-of-keys fix alone was not enough** - a union cannot declare a column no row carries - so the key is declared empty as well, and both the gaps-only and the mixed case are asserted |
| `mref` | 8 | Breadth dilution follows from having SEVERAL references, not from rank (D42). A species with two single-sequence assemblies is judged on the one its reads are actually on, and the substitution changes the TIER rather than only printing a note. A single-reference taxon gets no note at all. A DRAFT assembly blocks the substitution instead of letting one contig stand for a genome, and says the row may be under-called. With no refs manifest the question is left unasked |
| `stale` | 13 | Existence is not freshness (D45). A competitive bundle is rebuilt when its INPUTS changed, not merely when it is absent: a sidecar records the assemblies actually concatenated, keyed on accession, so a reference added, removed or swapped is caught with every mtime settled - the shape D42 hit, which mtime alone reports as current - while mtime still catches content changed under a known name. Missing and EMPTY are both stale. Indexes likewise: an index older than its bundle is the same disagreement one step on. Both decisive cases were run against the old existence-only predicate, which calls each of them current |
| `bundlecheck` | 13 | D45's predicate made ASKABLE (D61). It lived inside `build_pack_refs.py`, past the NCBI resolution and the subtree queries, and `--dry-run` returns before reaching it - so the one artifact stage 2 aligns against was the only database here with no read-only form. The check IMPORTS the predicate rather than copying it, and that is asserted on `__code__.co_filename` rather than by grep, because a second copy would answer with yesterday's rules while the builder rebuilt. It reads the bundles by default because only **10 of 90** carry the sidecar the predicate prefers - the other 80 predate it, leaving mtime, which a restore from backup defeats - so metadata alone answers *could not tell* for 89% of the pack, and `--quick` is asserted to say UNDETERMINED there rather than CURRENT. The content comparison is per taxid prefix against the manifest's `genome_size` and `n_contigs`, and the counting rule is tested FUNCTIONALLY against `fasta_stats` itself: a rule differing by the newline would be off a few bytes per contig and call every bundle on every deployment stale. The self-test puts a synthetic pack through 15 cases, including the two that metadata cannot see; three more mutations run against the REAL pack, where +1 bp in one manifest row exits 3. A missing neighbour and a de-duplicated one are reported with their consequence and do NOT exit 3, because both are permanent properties of this pack |
| `sketch` | 19 | D45's predicate arriving at the one artifact that had no check: the sylph sketch (D46). The shell `--check` is run against a STUB `sylph inspect`, so every verdict is exercised offline - a reference never sketched is caught with all mtimes settled (the state that actually shipped for 22 days), one newer than the sketch is caught by mtime, a sketched genome no longer in `refs/` is caught, and a database at the wrong `c` is stale too because it matches nothing while loading cleanly. Undetermined stays distinguishable from a verdict: an unreadable database exits 1, an absent one exits 3. The guard is asserted to be CALLED in the sylph branch, not merely defined; and no guard may count a list it deduplicated with Groovy's in-place `.unique()`, which made both this guard and the masking guard report "1 reference" for the twelve-file state they exist to report. A build RECORDS the genome set beside the database and `--stamp` records what the DATABASE holds rather than what `refs/` holds - derived from `refs/` it would agree by construction, the one thing it must be able to disagree with - and a record disagreeing with its database is caught, because a run trusts that file and cannot look inside the `.syldb`. The builder and the run are asserted to name the SAME record file: if those drift nothing fails, the builder writes a record nobody reads, and every run falls back to the timestamps-only check the record replaced. Finally the shipped artifacts: no reference newer than the sketch, its recorded genome count equal to the pack's reference count, and the record naming EXACTLY the pack's 375 references |
| `markerdb` | 23 | The same predicate at the marker database, whose input is `markers.tsv` (D47), and the mtime/set relationship is INVERTED from the sketch's: here the sets are blind to a marker **edited** under an id the manifest already knows - every set agrees while the sequence underneath is the one the last build retrieved, which is D41's wrong-gene shape - so mtime is the only signal that sees it, and the case is run with the mtime comparison removed, where it reports FRESH. A claimed tier the combined artifact does not hold is caught because that one fails UNSAFE: `marker_evidence.py` believes the manifest, searches, finds nothing and reports ABSENT, asserting the marker was looked for. Presence is compared and never sequence COUNTS - `concat_files()` de-duplicates within a marker while `n_nt_seqs` counts the per-marker archive, and the count-comparing version was run: it calls the correct shipped database stale on **116 of 166** markers, which is how a guard gets switched off. An index older than its source is caught, a missing one too; an absent database is stale (3) while an absent catalogue is undetermined (1), so a verdict stays distinguishable from a failure. `--check` is asserted to answer BEFORE `load_rankedlineage()`, because a read-only question costing minutes and gigabytes is one nobody asks twice. The guard is asserted to be CALLED in the marker branch, and the column parser to be a **top-level `def`**: written first as a closure in a local `def` and called as `idsOf(...)`, it failed to COMPILE under Nextflow 26.04's strict syntax - every run with markers enabled, which is the default - and that was found by extracting the function and running it, not by reading it. Finally the shipped database itself, checked in full on every commit, which is where the two build-integrity comparisons the launch guard leaves to `--check` are actually exercised |
| `envs` | 7 | The three build environments as FILES rather than a paragraph (they were documented only in prose, so an outside user had to reconstruct them from a description). These are curated specs, not `conda env export` output, and the invariant is why: the environment they were derived from still holds **kraken2 2.0.8-beta**, because installing kraken2 beside mmseqs2/hmmer/mafft makes the solver pick it - and it has no `--report-minimizer-data`, so an index built with it works while carrying none of the evidence stage 1 scores on. An export would therefore reproduce the trap. `btf-refs` is asserted NOT to declare kraken2 or bracken; both pins are asserted load-bearing rather than hygiene (kraken2 must match the version that BUILT the index or the format may not read back, and sylph's pin must equal the one in `modules/local/detect.nf` because a sketch must be readable by the version that queries it); and every file is asserted to declare only conda-forge and bioconda, since `defaults` imposes licence conditions on a reader and the originals had inherited unrelated channels from one developer's global config |
| `idxcover` | 6 | The Kraken2 index's question is COVERAGE, not freshness (D48): a watchlist taxon with no node in its pruned taxonomy can never appear in its report, and the arm still records `ok`. **19 of 90 are absent** on the shipped index - helminths, several protozoa, most dimorphic fungi - which follows from PlusPF's library selection and not from staleness, so a rebuild changes nothing. The number is pinned in BOTH directions: more is not automatically a regression and fewer is not automatically a fix, so either way the instruction is to re-measure and update D48 rather than adjust the number until the test passes. The bound is guarded harder than the gap: all 19 are sketched for sylph, which is on by default with `combining_rule: any_classifier_or`, so one arm is blind rather than the stage - and a taxon with NO node and NO sketched reference would be invisible to stage 1 entirely, which the check reports separately and this test fails on. Two overclaims are blocked: the pruned taxonomy (62,771 nodes) is never compared to the pinned taxdump's 3,016,439 - that comparison calls every correct index broken - and a present node is stated as necessary, not sufficient, since it can exist because a relative carried the genome. `taxo.k2d` is read with no kraken2 on the head node, and the layout arithmetic must add up to the byte or the reader REFUSES: a wrong stride does not fail, it reads the wrong 8 bytes of every record and reports them as taxids confidently |
| `extrarefs` | 11 | Extra references come from a MEASURED cover, not from curation rank (D44). Both input formats parse; a missing cover file is FATAL rather than silently empty, because loading nothing would build a one-reference-per-taxon pack while the operator believed a cover had been applied; the curated assembly stays first and keeps the plain filename earlier measurements used; pins follow in the cover's order and are not re-sorted, which would discard the measurement; and a pin already in the ranked list is not held twice - the two-rows-one-accession failure of D42, reachable again through this path, and verified to fail against the naive implementation |
| `refs` | 3 | Invariants over the shipped reference pack, all three violated on 2026-10-21 and each verified to fail against a reproduction of the broken state (D42): a `(taxid, accession)` pair appears ONCE; a taxon's row count equals its distinct accession count - `--max-per-taxon N` deduped candidates by FILENAME, so it re-downloaded accessions already held and credited a taxon with N references while holding one genome N times; and no competitive bundle is older than the genomes it is built from, because the FASTA and indexes are skipped whenever the output merely exists, leaving the manifest saying two references and the aligner using one |
| `matrix_blank` | 3 | A host-derived false positive is only suppressible if a control CONTAINS host (D14). *T. gondii* was `supported` in 18/18 lod_v1 samples off reads attribution showed were 100% human: the suppression machinery was never broken, the batch's only control was an NTC built at `host_fraction 0.0`, so `control_rpm` stayed 0 and `suppressed_control` could not fire. BOTH halves are asserted, because a test that checked only the matrix-blank case would also pass if everything were suppressed unconditionally |
| `neighbour_named` | 10 | `attributed_to_neighbour` must say WHICH organism took the reads (D18). It is the only tier that reports a watchlist organism absent on the strength of a different one, and it used to name that one by number - "reads favour neighbour taxid 620", where the analyst cannot see that 620 is *Shigella*, let alone that the pack's own notes say the pair is inseparable at species level. The fallback is asserted too: the name comes from a manifest the scorer must not depend on, so without it the bare taxid returns and scoring is unchanged |
| `contest` | 8 | A neighbour may only take a finding the target has no claim on. The attribution branch never looked at the target's own evidence, so a target covering 0.667 of its reference lost to a sibling at ratio 0.9032 and was reported absent - and since that sibling is not a watchlist taxon, NO finding anywhere carried its 3.3 M reads. THE NOT-CONTESTED CASE IS WHAT KEEPS THIS HONEST: D18 measured 25 false positives from applying a neighbour gate to every attributed row, so a weak target's attribution is left exactly as it was |
| `reachable` | 6 | A perfect finding at each plausible read length must reach the top tier - **the test the suite did not have, and its absence cost the project its whole top tier on real data**. `min_alignment_length` was a flat 100 while a 2×101 library trims to a mean of 98-99, so 81 of 82 rows in the Zymo isolate batch failed it and every `confirmed_in_silico` call in nine runs was on 150 bp SIMULATED reads. The suite could not see it because its one `mean_alignment_length` fixture is exactly 100.0, the single value that passes. Every other test asserts a gate REJECTS bad evidence; this one asserts gates ACCEPT good evidence, parameterised over the axis the benchmark holds fixed |
| `ctrlaudit` | 14 | Control subtraction must be readable off a table, not regexed out of prose (D21). Across six batches it moved 108 rows with the `reasons` column as the only record, hiding three things: which control decided a suppression, that 11 suppressed rows carried MORE of the organism than the control and were held down purely by the 5× margin, and that a control can hold a select agent at a level that suppresses nothing today while standing as a floor the moment a sample lands just under it - lod_v3_threearm's NTC holds *B. anthracis* at 32.58 RPM and nothing said so. Each world gets its own coverage file, because mixing them poisons run-global `control_rpm` |
| `crossbleed` | 4 | D16a's mechanism, kept as a test rather than as a batch on a disk. `simulate_benchmark.py` now refuses to build a batch mixing a hopped control with other scenarios, so the demonstrating batch is unconstructible - which makes this fixture the only standing record of what is at stake: a control from one scenario suppressing genuine detections in another, because `control_rpm` is keyed on taxid alone across the whole run. The numbers are the measured ones: a 1,000 RPM control row turns a 1,816.7 RPM `supported` call into `suppressed_control` |
| `hopguard` | 5 | A deliberately contaminated control cannot share a batch (D16a). The scenario is correct alone and incoherent in company - its NTC holds real target reads by design. Refused at the batch level rather than scoped in the scorer, because that run-global reach is exactly what D14's matrix blank depends on: re-scoring lod_v4 without its controls turns 15 *T. gondii* and 15 *E. coli* rows back into false positives. Asserted in BOTH directions, since a guard that refuses everything passes the first check and is the likelier mistake |
| `protfloor` | 5 | The protein tier must use the floor the curator wrote, not a blanket discount (D23). `aa_min_id = max(0.50, min_id - 0.15)` subtracted 15 identity points from every marker regardless of curation, overwriting exactly the distinction the catalogue encodes: housekeeping genes that cross-match across genera carry 0.95-0.98, and the divergent-virus block carries 0.78-0.90 ON PURPOSE with a note that the HMM tier is what carries detection there. Both directions from one catalogue, because a change that merely tightened everything would pass the first assertion - and tightening the filovirus markers is the specific harm |
| `markerexp` | 10 | The second expectation axis, end to end through the CLI: `--marker-states`, `load_marker_states` reading a real `marker_evidence.py` file, and the exit code that gates CI. The axis exists because a findings row is only ever written for a WATCHLIST taxon - divergent_clostridium_bont spikes *C. butyricum*, deliberately a near neighbour, so its only expectation could be a `max_tier` upper bound that "no finding" satisfies. The scenario was green on every run while the toxin it exists to test went unreported |
| `tiergap` | 11 | A tier that was never built must not be reported as having been searched. `tiers_supplied` answers "did the run get an input file for this tier" and was used as if it answered "does this marker have anything in it", with `tried` appended unconditionally - so a marker whose HMM was never built emitted `tier_used = dna,hmm,protein` and "searched via dna,hmm,protein; no tier met thresholds". On the shipped pack that was 13 markers, almost always the HMM tier and every one a divergence-catching marker on a high-consequence agent; orthopox_ATI and mpox_C3L were missing the NUCLEOTIDE tier instead, which is worse - a Variola marker claiming a nucleotide search that never happened |
| `markerjoin` | 4 | Three hand-maintained files agree on literal strings, so something must check: `markers.tsv.taxon`, `marker_rules.yaml`'s keys and `watchlist.resolved.tsv.name` name the same organisms. Drift degrades SILENTLY - a rule block whose key stops matching simply stops gating, which is indistinguishable from a taxon that has no rules, which is D17a's permissive path |
| `leadcap` | 10 | A family cannot be confirmed as an organism, however good the geometry. A watchlist entry resolved above species exists so a divergent member with no reference of its own still surfaces; what that supports is "a member of this group is here", which is a lead. The cap is explicit rather than left to breadth arithmetic, because a family bundle's breadth divides the covered span by the summed `ref_bp` of every genome sharing the taxid - so whether it clears `min_occupancy` depends on how many genomes the family happens to have, and adding one would change the verdict |
| `leadjoin` | 6 | `resolve_watchlist` writes the descendants file and `evaluate_benchmark` reads it: two files, one hand-written column contract, and a rename on either side degrades SILENTLY - `load_covers()` returns an empty mapping, every group lead reads as an unlisted false positive, every species it covered goes back to a flat miss, with no error anywhere. The evaluator's self-test checks the join logic against a hand-built dict; this checks that a real file feeds it |
| `defaults` | 5 | An argparse default that disagrees with `nextflow.config` is a trap for re-scoring. D18 set `min_neighbour_uniformity = 0` in the config after measuring the gate on real reads; the argparse default stayed 0.50, so every manual re-score that forgot the flag silently reproduced PRE-D18 behaviour - which happened during the D16a review and produced 7 tier changes where the shipped configuration produces 5, reading as an error in `decisions.md` that was not there. Only thresholds present on both sides are compared, and the last check fails if the pairs stop matching the files at all, so this cannot degrade into testing nothing |
| `sketchauth` | 5 | The sketch's authority is the manifest, not the directory it reads (D59). `build_sylph_db.sh` sketched whatever FASTA was in `refs/` and `--check` compared the sketch against that same directory, so both sides of the comparison were the same wrong thing. Measured while reverting a reference set: 31 files the manifest no longer named stayed on disk, the rebuild sketched them, the sketch reported 406 genomes against a 375-reference manifest, and `--check` answered FRESH - every count agreeing with itself |
| `metabuli` | 9 | The last unchecked arm; opt-in is not the same as covered (D60). The Kraken2 test does NOT transfer: that one asks whether a watchlist taxon has a node in the index's PRUNED taxonomy, discriminating because `kraken2-build` keeps only nodes spanning sequence it loaded, whereas a Metabuli database copies the WHOLE taxdump - so every taxon has a node and the same test would pass unconditionally. **A check that cannot fail is not a check.** The discriminating input is `taxID_list`, the taxids the database holds sequence for, asked over the clade |
| `roots` | 15 | No script may name a path on the machine it was written on (D51). Scanned: EVERY tracked file, comments and documentation included, for an absolute path under a deployment volume or a home directory, with `docs/decisions.md` the one exemption because its entries cite the paths they were measured on. Two earlier versions were weaker and both were defeated - the first anchored on `=`, `default=` and `:-`, which a default reached through a conditional expression evades; the second exempted comments. Plus the loop that runs every builder and read-only check with NO environment and requires each to refuse by name, because a `--check` reading a root that is not yours reports a healthy pack as having no database at all |
| `lock` | 14 | The pack is reproducible only if the assemblies are PINNED, not re-chosen (D52). The rule is deterministic; NCBI's catalogue is not, so two people build two different packs from one watchlist with no error on either side. The checksum is of the DECOMPRESSED sequence, which is load-bearing: the container md5 reported a replaced assembly for byte-identical sequence fetched by another route, and a guard that cries wolf at every user gets switched off. Also that the lock OUTRANKS the disk (without it a built pack satisfies `--max-per-taxon` from its own files and a corrected lock is inert), that an unsatisfiable lock FAILS the build rather than falling back to latest, and that locked-absent stays distinguishable from not-in-the-lock. Offline, so the four live NCBI behaviours are greps protecting paths verified by hand on 2026-10-23 |
| `install` | 37 | A clone has to be able to get the databases, and to know what it got (D53). The installer is the only part of this project every user runs before anything works. The manifest carries every column the installer reads and names only routes it knows; only the artifacts that CANNOT be re-created are published; the checksum is verified BEFORE anything is unpacked, because a truncated taxdump does not crash `resolve_watchlist.py` - it resolves fewer taxa and reports success; a `PENDING` url is refused rather than guessed at; the build order (refs → host → masking → sketch) is pinned by a test because both orderings were learned from a failure; and the disclaimer travels with the data, since an archive that arrives without it reads as a curated reference resource |
| `container` | 30 | The container profiles declared no image at all until v0.4.0 (D54). `-profile docker` and `-profile singularity` sat in `nextflow.config` for months while not one process declared a `container`: the profile looked supported and was absent. `conf/containers.config` and `envs/process/*.yml` are GENERATED from each process's own `conda` directive, because the hand-written version is the one that rots - a pin moves, nobody updates the map, and the two routes run different versions of the same tool, correct one way and wrong the other, visible only in results |
| `wiki` | 11 | The wiki is maintained in the repository, so it can be checked like code (D56). A 63 KB README is a reference manual on the front page: nobody reads it to the end and it goes stale invisibly because no test touches prose. The split only helps if the pages and their links exist - a sidebar entry pointing at a page nobody wrote renders as a live link to an empty page, the documentation version of a profile that declares no container. Every page the sidebar and README promise, every wiki link, every repository link as an absolute URL, every image pointing at a file that exists, and the workflow that syncs deletions too |

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
  `build_pack_refs.py` as a post-write step. ~~**Not automated:**~~ **Now automated
  (D61)**, and no longer only on rebuild: `scripts/check_pack_bundles.py` re-reads every
  bundle and re-applies that post-write check plus the unprefixed-header one, so the 80
  bundles built before either existed have now been subjected to both. Neither fires, and
  all 90 bundles and 180 indexes match the manifest per taxid prefix (2026-10-26, 73 s).
- **Taxonomy drift**: the previously-default snapshot measured 201 days old, +288,764
  taxids behind current, 2,581 descendant nodes short inside watchlist clades.
- **Masking**: 20,040,354 bp over 375 references (0.348%); 104.2 kb vector across 45
  contigs. As of 2026-10-21 all **375** pack references that have sequence carry a result;
  the twelve Filoviridae group-entry files had never been masked (D34, when the pack held
  371) and came back at 0 bp. The script now
  reconciles the whole pack and exits non-zero when a reference has no result, so the
  count no longer has to be checked by hand. **And in the other direction, from 2026-10-26
  (D62):** 31 BEDs and 31 `mask_manifest.tsv` rows belonged to no reference the manifest
  named — D58's reverted cover, left behind because the manifest merge carries every old
  row forward — so the artifact recorded 20,148,007 bp over 406 results while this figure,
  counted over the manifest, was right. Pruned; `mask_manifest.tsv` now holds exactly 375
  rows summing to 20,040,354 bp. `mask_pack_refs.py` reports orphans on every exit path and
  exits non-zero only when one names a contig still in the pack, since `coverage_stats.py
  --mask-bed` reads the whole directory into one contig map and such an interval masks live
  sequence; none of the 31 did.
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
