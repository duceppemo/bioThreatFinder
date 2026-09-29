# Tool selection and design decisions

One entry per decision: what was chosen, the evidence, and **what would reverse it**.
Recorded because the reasoning is the reusable part — a tool list without its evidence
rots in a year and nobody knows whether it was measured or guessed.

Evidence is tagged: `[IND]` independent peer-reviewed · `[TOOL]` the tool's own paper ·
`[SRC]` read from source code · `[MEASURED]` measured in this project · `[REASONED]` a
judgement call with no measurement behind it.

**Nothing in this project has been benchmarked against ground truth yet.** Where a
decision rests on `[REASONED]`, treat it as provisional.

---

## D1 — Two stages: permissive screen, strict confirmation

**Chosen:** stage 1 raises candidates permissively; stage 2 decides by alignment.

Stage-1 classifiers are built for community profiling, not for "is organism X present".
A handful of k-mers shared with an environmental relative puts a high-consequence taxon
in a report. `[IND]` Wright et al. 2026 measured Kraken2 unfiltered at **264 true vs
3,938 false** taxa per sample (*Microb Genom*, doi:10.1099/mgen.0.001739).

Supported from an unexpected direction: `[IND]` WOAH Terrestrial Manual **Ch. 1.1.7**
states that *"confirmatory adjunct or secondary diagnostic assays need to be validated
only for their analytical performance … and not to the full diagnostic extent"* — a
standards body formally sanctioning screen-then-confirm.

**Would reverse it:** evidence that stage-2 confirmation cannot be made sensitive enough
at trace abundance, making the second stage a pure sensitivity cost.

---

## D2 — Kraken2 primary, KrakenUniq opt-in ⚠️ *reversal of an earlier decision*

**Chosen:** Kraken2 with `--report-minimizer-data` as the primary nucleotide arm.
KrakenUniq available but off by default.

**This reverses my own earlier recommendation**, which was wrong. I had argued KrakenUniq
should be primary because "Kraken2 only counts reads". That is untrue for ≥2.1.0, and the
reason is stronger than a feature overlap:

- `[SRC]` Kraken2's `src/readcounts.h` carries the header *"originally developed for the
  KrakenUniq taxonomic classification system"* — it is literally KrakenUniq's code, by
  the same author, using the same `hyperloglogplus.h` at the same default precision
  **p=12**. The only difference is counting distinct **minimizers** instead of distinct
  **k-mers**: a resolution difference, not an algorithmic one.
- `[IND]` Oskolkov 2026 (*Front Microbiol*, doi:10.3389/fmicb.2026.1603339) is the only
  published comparison and finds **comparable F1** once the threshold shifts by the 2–3×
  that minimizer subsampling implies: ~1,000 unique k-mers (KrakenUniq) vs **~2,000–3,000
  distinct minimizers** (Kraken2). It concludes the k-mer filter alone is the most
  effective FP filter, better than read count or coverage.
- KrakenUniq has had **no commit since December 2023**; Kraken2 shipped 2.17.2 on
  2026-09-23. For a pipeline to be maintained for years that asymmetry outweighs
  "slightly higher accuracy" at one calibration point in one domain (that paper's domain
  is ancient DNA, not biosurveillance).

Cost avoided: ~1400 GiB transient disk, a 250–450 GiB resident index, a **Jellyfish 1.x**
build dependency this machine cannot satisfy, and 12+ hours.

Three caveats on Kraken2's columns, `[SRC]`: columns 4–5 are **clade-aggregated**, so
threshold at species rows not genus rows; minimizers are keyed by the minimizer's
database LCA rather than the read's final call; and upstream issue #445 (column 5
exceeding `kraken2-inspect`) is still open, so treat column 5 as a monotone score to
calibrate, not a literal breadth measure.

**Would reverse it:** the bake-off (`scripts/run_bakeoff.sh`, variant
`plus_krakenuniq`) showing a materially better confirm-PPV that justifies the cost. A
prebuilt KrakenUniq index exists, so even then **do not build one**.

---

## D3 — Metabuli kept despite poor independent precision

**Chosen:** keep Metabuli as a **recall** arm, never as a primary caller, pinned to a
specific version and database snapshot.

Two reviews appear to conflict; they do not, once its role is stated correctly.

- Against: `[IND]` LEMMIv2 (*Genome Biology*, April 2026, doi:10.1186/s13059-026-04089-9)
  ranks Metabuli 1.0.2 **lowest on precision** alongside Centrifuger. `[IND]` The CAMI
  portal (*NAR* 2025) ranks Metabuli 1.0.8 **second at species level** on short reads.
  The disagreement tracks *whether the benchmark controls the reference database and
  calibrates thresholds per tool* — which matters more for our design than either
  ranking.
- For: `[IND]` Perseus (*Bioinformatics* 2026, doi:10.1093/bioinformatics/btag687)
  measured Kraken2's false-assignment rate at **~50% in family-level exclusion** on long
  reads. Divergent reads are not left unclassified — they are **confidently assigned to
  the wrong known taxon**, frequently on a median of five k-mers. For surveillance that
  is the worst available error shape: the novel agent is not missed, it is mislabelled as
  something familiar and dismissed. Metabuli's joint DNA+amino-acid k-mers are a
  documented mitigation.

Precision is stage 2's job. Recall at stage 1 is what cannot be recovered later.

**Known operational risks** (open upstream issues, checked directly): broken database
downloads (#195, #190, #185), taxid inconsistency in reports (#193), **segfault on a
mismatched taxdump** (#175), and a trailing tab on unclassified rows that breaks TSV
parsing (#194 — handled in `parse_screen.py`).

**Would reverse it:** the bake-off showing `plus_metabuli` adds nothing on the
`divergent_*` scenarios. A lower-risk route to the same capability is Kaiju v1.10.3 or
DIAMOND-based translated search.

---

## D4 — Any-classifier OR, not consensus

**Chosen:** any arm raises a candidate. Concordance is **recorded**, never required.

`[IND]` The one rigorous 2025–2026 study on combining classifiers (*PLOS Comp Biol* 2026,
doi:10.1371/journal.pcbi.1013883) found **22% species concordance** between Kraken2/Bracken
and MetaPhlAn4 — and 61% of Kraken-unique species *were* in MetaPhlAn's database yet
undetected, so disagreement is not merely a database artifact. That figure argues a union
rule would be a false-positive problem and an intersection rule would be dominated by
MetaPhlAn's near-zero recall in low biomass. We take the union at stage 1 precisely
because stage 2 supplies the precision.

That study optimises for **discovery recall in community profiling**, which is not our
objective — noted so it is not over-read. **No 2025–2026 evidence supports
multi-classifier consensus for precision in pathogen detection.**

**Would reverse it:** measured evidence that a specific combining rule beats
OR-then-confirm on confirm-PPV at equal sensitivity.

---

## D5 — Reference masking before threshold tuning

**Chosen:** soft-mask vector, host, low-complexity and cross-kingdom regions; exclude
masked bases from breadth denominators; report reads landing in them.

Two independent reviews put reference curation **ahead of assembly** as a specificity
lever: `[IND]` 40–94% fewer false detections with 100% true-positive preservation
(*J Clin Microbiol* 2026). `[IND]` Wright et al. traced an *H. pylori* call with **zero**
mapping reads to **phage sequence inside the reference genome**.

Soft- not hard-masking is the load-bearing choice: hard-masking makes that failure
*invisible*, because the reads simply stop aligning and the taxon vanishes with no record.

`[MEASURED]` on this pack: 20,028,976 bp masked (0.349%) — 93.9 kb vector across **44
contigs**, 79.7 kb cross-kingdom across 49, 2.7 kb host across 5. Far below the headline
figure, correctly so on a 281-of-363 reference-grade set, but 44 contigs carrying vector
sequence is the exact Wright mechanism.

**Would reverse it:** nothing plausible. It is cheap and its failure mode is visible.

---

## D6 — Enrichment is NOT the biggest lever here ⚠️ *reversal of an earlier claim*

**Chosen:** shotgun only. If enrichment is ever added, **always retain an unenriched arm.**

I previously told the user enrichment was the single biggest sensitivity lever, "orders of
magnitude" ahead of software. **For a mission that must catch divergent agents that is
wrong.** `[IND]` Yamaguchi et al. 2025 measured:

- capture needs **>65% nucleotide identity over ≥1 kb** to enrich at all;
- **off-panel agents get worse**: simian haemorrhagic fever virus read abundance fell
  **25-fold** after capture; Nodamura virus segments 33- and 50-fold;
- panels age — 2017–18-designed influenza probes gave 40–86% coverage for H5N1 clade
  2.3.4.4b vs >98.5% for a 2024 panel.

Enrichment is a large lever for **on-panel** targets in low-biomass matrices and an
actively **negative** one for the novel agent this tool exists to catch. The published
answer if both are wanted is **teamNGS**: pool enriched and unenriched libraries 90:10.

Adaptive sampling is likewise 1.2–5× realised for 50–80% yield cost — free polishing, not
a strategy.

**Would reverse it:** a bounded, frequently-refreshed panel and an explicit decision that
on-panel sensitivity matters more than novel-agent detection. That is a mission change.

---

## D7 — Watchlists and markers keyed by name, not identifier

**Chosen:** names resolved against a pinned taxdump at build time; unresolvable names are
a hard error.

`[MEASURED]` It paid for itself immediately: 89/89 resolved while flagging 18 changed
virus ranks, 11 ICTV renames resolving via equivalent names, and 4 malformed neighbour
entries. Hardcoded taxids would have rotted silently, and a rotted taxid means an
organism stops being screened for **with no error**.

`[MEASURED]` The same principle had to be applied twice more after being missed:
`rankedlineage.dmp` carries only current scientific names, so family-level relaxation
found nothing for "Marburg marburgvirus" (now a synonym) until the lookup was made
synonym-aware.

**Would reverse it:** nothing. This is the highest-value engineering decision in the
project.

---

## D8 — Pinned, dated taxonomy snapshot

**Chosen:** a project-owned dated snapshot; the shared `/media/36tb/db/taxdump` left
untouched.

`[MEASURED]` The shared snapshot was **201 days old**. Drift to current: **+288,764
taxids** (+10.5%), 2,906 merged/removed. Two silent consequences: `kraken2-build`
downloads *current* library sequences whose taxids would not exist in an old taxonomy, and
watchlist clade coverage was 2,581 descendant nodes short (+1,969 *Salmonella enterica*,
+253 *L. monocytogenes*, +150 Influenza A, +36 *F. tularensis*). All 89 taxa still
resolved to the same taxids, which is precisely why it was invisible.

**Would reverse it:** nothing. Guards now refuse a snapshot >90 days old.

---

## D10 — Read-identity floors are per platform, not only per group

**Decision:** `min_read_identity` resolves platform+group → platform → group → default,
and the findings row names which rule was applied. Illumina keeps 0.98 for
bacteria/fungi; ONT starts at 0.90 (virus 0.80, parasite 0.88). The marker tier gets the
same treatment via `--platform` / `--ont-identity-slack`.

**Why:** `[MEASURED]` `bin/score_evidence.py` contained no platform term at all. A 0.98
read-identity floor is a statement about **sequencing error** as much as about
phylogenetic distance. Illumina substitution error is ~0.1–1%, so 0.98 leaves room to
discriminate species. ONT does not reach it: R10.4.1 simplex means ~97–98% with current
basecallers, and `badread`'s default — what our own ONT benchmark simulates — is a 95%
mean, verified from `badread simulate --help`. So **every bacterial and fungal ONT finding
failed the identity gate on sequencing error alone**, at every spike level, and capped at
`supported` with a reason string indistinguishable from genuinely weak evidence. The
`plumbing_positive_control`, whose stated purpose is to fail loudly on misconfiguration,
would have failed permanently on its ONT arm for this reason.

**What this costs:** on ONT, read identity can no longer discriminate near neighbours.
That burden moves to competitive mapping (`attributed_to_neighbour`) and the marker tier,
and the identity gate degrades to a genus-level sanity check. That is a genuinely weaker
guarantee on that platform, and it is stated in `pack.yaml` rather than hidden in a number.

**Status: PROVISIONAL.** The ONT numbers are reasoned from platform error rates, not
measured on this pipeline. They are the first thing to set once the benchmark runs.

**Would reverse it:** measurement showing ONT near-neighbour confusion at 0.90 that
competitive mapping and markers do not catch. Then the floor rises and ONT confirmation
depends more heavily on the marker tier.

**How it hid:** two regression tests asserted the floor was **strict**
(`>= 0.98`) rather than **reachable**. They now assert reachability per platform, and a
fixture-level test scores the same evidence on both platforms.

---

## D11 — The two stage-2 gates are mutually incoherent, and that is recorded rather than patched

**Finding:** `[MEASURED]` `min_unique_reads` and `min_breadth` are ANDed, and breadth is
bounded by read count — with reads of length *L* on a genome of size *G*, *N* reads can
reach at best *N·L/G* breadth. So whenever `min_breadth · G / L > min_unique_reads` the
read floor **cannot** bind: the stated limit of detection understates the real one by
exactly that ratio.

Measured on the shipped human pack at 150 bp, using genome sizes from the built refs
manifest: **75 of 89 taxa** have a decorative read floor. Worst cases —
*Schistosoma mansoni* 243×, *Angiostrongylus cantonensis* 196×,
*Echinococcus multilocularis* 91×. *Balamuthia mandrillaris*, whose own watchlist note
says "low read counts are expected and still meaningful", needs ~3,159 reads to clear its
own breadth floor against the 50 the file declares. In the other direction breadth never
binds for *R. prowazekii* or Influenza A, so neither column does the job it claims.

**Decision: report, do not recalibrate.** Recomputing one gate from the other requires
choosing a spread factor, and choosing it from reasoning rather than measurement would
replace one arbitrary set of numbers with another — the exact thing this project is trying
not to do, and these thresholds are already slated to come from the benchmark.
`scripts/check_threshold_coherence.py` prints the per-taxon table and
`tests/test_pipeline.py` **pins the counts at 75 of 89**, so a pack edit shows up as
coherence improving or regressing instead of silently staying broken.

**Would resolve it:** the first real benchmark run. Per-taxon detection rates against
declared rpm give the actual LoD, and both gates get set from that.

---

## D9 — Rejected for now

| Option | Why not yet |
|---|---|
| **Assembly tier** | `[IND]` Recommended by the mLife 2026 clinical review and it raises contig-level specificity, but reviews put it *behind* reference curation, and read-level calling must not be gated on assembly succeeding. Roadmap. |
| **Empirical null / calibrated FDR** | The right end state. Needs a background panel of negative/unrelated metagenomes we do not have. Fixed thresholds are the honest interim. |
| **Protein language models** (LucaProt, Rider, LucaVirus) | `[TOOL]` Competitive on RNA-virus recall, but **every published figure is the developers' own** and no independent head-to-head exists. A 2026 preprint also shows pLM embedding space models viral proteins *worse* than its general benchmarks imply. Not a decision-maker. |
| **Kun-peng, Perseus** | `[IND]` Both peer-reviewed 2026 and both directly relevant — Kun-peng claims fewer FPs than Kraken2/Centrifuger/KrakenUniq in 4.1 GB RAM; Perseus abstains or backs off a rank rather than over-assigning, which matches our tier philosophy. Watch list. |
| **Multi-reference / pangenome bundles** | Single reference understates a divergent strain. Correct fix, not yet built. |

---

## Standing principle, learned the hard way

**Well over thirty** distinct bugs in this project have shared one shape: *a mismatch that
produces an empty or degraded result indistinguishable from a clean negative.* Rotted
taxid · GTDB namespace · refs base path · crossgroup subset · duplicated bundle sequence ·
stale taxonomy · kraken2 2.0.8 lacking the minimizer flag · marker manifest truncation ·
refs manifest truncation · sylph filename keys · ART fold coverage · a silent Kraken2
death · the ONT identity floor (D10) · two genomes never masked but recorded as clean ·
minimap2's ignored exit code · a GTDB index stamped `NCBI` because the stamp read a flag ·
`--with-maps` overwriting the taxonomy it had just pinned · seven threshold params
declared and wired to nothing · a header-only samplesheet exiting 0 · an index-hop that
moved 0.02 target reads · a benchmark axis declared in YAML and read by no code · a
bake-off variant identical to its own baseline · a vaccine check that could never fire.

So the first question for any new component is **"could this silently produce nothing?"**
— and the answer must be a guard that fails loudly, plus a check against the written
artifact rather than an assertion in code.

Three sharper lessons, each paid for:

1. **A guard that asserts how a script was invoked is not a guard.** The GTDB namespace
   stamp, the Bracken read-length model, the KrakenUniq HLL precision and the host-genome
   path check all passed their own logic while failing their stated purpose. Check the
   artifact.
2. **A threshold must be reachable, not merely strict.** Tests that assert a floor is
   high lock in a floor that nothing can meet (D10, D11).
3. **Stamps make a resume safe sequentially and say nothing about concurrency.** Two
   builders sharing a database directory can interleave partial downloads, and a corrupt
   reference in an index does not error — it misclassifies. This was learned by doing it.

Several fixes were themselves vacuous on the first attempt and were caught only by
re-reading the written artifact.
