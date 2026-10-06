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

## D12 — `neighbour_ratio`'s denominator is 98% background cross-mapping and 0% neighbour ⚠️ *confirmed by read attribution*

**Finding:** `[MEASURED]` 2026-10-03. The gate is documented as "do these reads favour the
target or its near neighbour?" (`coverage_stats.py:39`). By direct read attribution it is
not measuring neighbour competition at all.

The simulated reads carry their source accession in the read name, so every alignment in
the competitive bundle can be traced to the genome it actually came from. For
`lod_gram_negative`, taxid 632, MAPQ>=10, reads assigned to the NEIGHBOUR reference
(taxid 630, *Y. enterocolitica*):

| sample | neighbour reads | from *E. coli* | from *Y. pestis* | % *E. coli* |
|---|---|---|---|---|
| r1 rpm10 | 2,763 | 2,700 | **0** | 97.7% |
| r1 rpm100 | 2,731 | 2,664 | **0** | 97.5% |
| r1 rpm1000 | 2,708 | 2,660 | **0** | 98.2% |
| r2 rpm100 | 2,773 | 2,718 | **0** | 98.0% |
| r3 rpm100 | 2,822 | 2,749 | **0** | 97.4% |
| r1 rpm10000 | 2,737 | 2,673 | **0** | 97.7% |

**Not one read from the spiked target is assigned to the neighbour, in any sample.** The
remaining ~2% is *H. haemolyticus* and *N. lactamica*, also background. The reads are
*E. coli* O157:H7 Sakai (NC_002695.2 plus its plasmids), which is the standing background
at 2%. There is no Yersinia neighbour in these samples at all, so the "competition" the
denominator represents does not exist.

**The numerator has the same defect, and at low abundance it is worse.** Reads assigned to
the TARGET reference, by true origin:

| level | target reads | true *Y. pestis* | *E. coli* | % false |
|---|---|---|---|---|
| 10 rpm | 427 | 225 | 202 | **47.3%** |
| 100 rpm | 2,336 | 2,114 | 222 | 9.5% |
| 1000 rpm | 22,320 | 22,128 | 191 | 0.9% |
| 10000 rpm | 221,439 | 221,209 | 230 | 0.1% |

So at 10 rpm **nearly half of what is counted as the target is *E. coli***. That corrupts
breadth, `depth_cv` and `mean_identity` too, not only the ratio.

**Both pedestals are constant, and their asymmetry is the whole bias.** *E. coli*
contributes ~2,700 reads to the neighbour reference and ~210 to the target reference, at
every spike level. The neighbour recruits **~13x more** background than the target does -
a property of the two REFERENCES, not of the sample. So in this background

    neighbour_ratio  ~=  (true_target_reads + 210) / 2,700

and clearing the 1.5 threshold needs about **3,840 true target reads** no matter what, with
or without a neighbour present. `min_neighbour_ratio` is therefore an undeclared read floor
whose height is set by background composition and by which reference happens to be more
*E. coli*-like.

**Falsifiable prediction, with the interval the replicates actually support.** Reads scale
linearly with spike level, so 3,840 true reads is reached near 180 rpm - which is why the
measured LoD came out as 1000: the series steps 100 -> 1000 and brackets the real threshold.
Computed per replicate from `-F 0x904` data rather than from one anchor (true reads at
100 rpm 2,114 / 2,233 / 2,166; neighbour pedestals 2,725 / 2,694 / 2,710; target pedestals
206 / 220 / 207):

    r1  (1.5 x 2725 - 206) / 21.14 = 183.6 rpm
    r2  (1.5 x 2694 - 220) / 22.33 = 171.2 rpm
    r3  (1.5 x 2710 - 207) / 21.66 = 178.1 rpm

So the prediction is a **band of 171-184 rpm, centre ~178**, not a point. An earlier version
read "~180 rpm... if it is found anywhere else, this arithmetic is wrong", which would score
a cliff at 170 or 190 as a refutation of arithmetic that in fact predicted it. At the cliff
the per-sample neighbour spread straddles the threshold - (3840+210)/2811 = 1.44 fails while
(3840+210)/2581 = 1.57 passes - so with three replicates the cliff is a band by construction.

### CONFIRMED 2026-10-06 — the cliff is at 176.6 rpm, inside the predicted 171-184 band

`runs/d12fine`: 120/150/170/190/210/250 rpm x3 replicates, ART at 150 bp / insert 350, batch
gated on 183 ART commands and **0 iss** before the pipeline was allowed to start (condition 1).
Pipeline exit 0.

| level | tier | mean neighbour_ratio |
|---|---|---|
| 120 | 3/3 `attributed_to_neighbour` | 1.038 |
| 150 | 3/3 `attributed_to_neighbour` | 1.290 |
| **170** | **3/3 `attributed_to_neighbour`** | **1.448** |
| **190** | **0/3 - all `supported`** | **1.604** |
| 210 | 0/3 | 1.722 |
| 250 | 0/3 | 2.087 |

Linear interpolation of the 1.5 crossing between the bracketing levels gives **176.6 rpm**,
against a prediction of 171-184 with centre ~178 - inside the band and 1.4 rpm from its
centre. The per-replicate spread behaves as predicted too: at 170 rpm one replicate sits at
1.4989, which is the threshold to four decimal places.

**Condition 3 is satisfied explicitly, not by inference.** The `reasons` strings name the
deciding gate at each side of the cliff:

    170 rpm: reads favour neighbour taxid 630 (ratio 1.4989 < 1.5)
    190 rpm: read threshold met; failed: depth_cv 3.62 <= 2.50

So `neighbour_ratio` is what binds below the cliff and releases above it, and `depth_cv`
takes over afterwards - which is why the tier above 190 rpm is `supported` rather than
`confirmed_in_silico`. That is consistent with the prediction, which was only ever about the
neighbour gate's crossing.

**What this establishes.** The mechanism is no longer inferred from one batch's pedestals: a
quantitative prediction made before the data existed, from read attribution on a different
batch, located a threshold crossing to within 1.4 rpm. `min_neighbour_ratio` is therefore an
undeclared read floor whose height is set by background composition - about 3,840 true target
reads in this background - and the measured "limit of detection" of 1000 rpm was an artefact
of a decade grid straddling it. The real detection limit is far lower than the LoD series
could show.

**Not affected by the 2026-10-06 `samtools depth -aa` correction.** That bug inflates breadth
and occupancy on FRAGMENTED references; *Y. pestis* has 4 complete replicons and every one was
emitted under `-a`, verified by re-running the patched `coverage_stats.py` over these bundles
(breadth and neighbour_ratio identical to the last decimal, depth_cv differing only in the
fourth). The numbers above stand.

**Three conditions the test batch must meet, or the result cannot be interpreted:**

1. **It must be ART-simulated.** `simulate_benchmark.py:2416` defaults `--illumina-simulator`
   to `auto`, and `auto` PREFERS InSilicoSeq when installed - which it is, in `btf-sim`.
   `lod_v1` was built with an explicit `--illumina-simulator art` (138 ART commands, zero
   iss). iss has NO insert-size parameter and a model-fixed read length, and insert alone
   moves these coefficients 10-28%, so a cliff measured on an iss batch cannot be compared
   against this ART-calibrated prediction. Running the finer grid with default flags would
   therefore refute the arithmetic for a reason that has nothing to do with it.
   `check_pedestal_model.py` now detects an iss batch and says exactly this.
2. Same read length (150 bp) and insert (350) as the coefficients were measured at.
3. **The cliff must be attributable to `neighbour_ratio`.** Detection requires EVERY gate,
   and D14 documents a case where `depth_cv` binds alone. If another gate binds first near
   180 rpm, the cliff moves for reasons this arithmetic says nothing about - so the
   per-sample `reasons` string must be read before the result is called a refutation.

**Masking does not currently touch this.** `mask_pack_refs.py` has a `crossgroup` source,
but it compares "only across distant groups" (bacteria vs virus/fungus/parasite) by design,
so bacteria-vs-bacteria recruitment is invisible to it: `masked_bp` crossgroup is **0 bp**
for both 630 and 632. Masking conspecific-family recruitment is therefore an available
lever that has never been pulled, and unlike a threshold change it would fix the numerator
and denominator together.

**Decision: still do not patch the metric; the fix is now a choice between three named
options, and the measurement to choose between them is cheap.** What changed since this
entry was first written is that the MECHANISM is no longer in doubt - it is background
recruitment, quantified, with zero neighbour contribution. What is still undecided is the
remedy, and the options differ in what they would cost every one of the 89 taxa:

1. **Mask within-family recruitment.** Extend `crossgroup` to same-group genomes, or add a
   `background` source masking against the pack's own standing background. Fixes numerator
   and denominator at once and needs no threshold. Risk: over-masking a genuine shared
   locus, which costs sensitivity silently - the thing masking is supposed to prevent.
2. **Require the neighbour to be credible before it can out-vote.** Gate attribution on the
   neighbour's own breadth and occupancy, not its raw count. A 2,700-read pedestal spread
   thinly over 4.5 Mb has negligible breadth, so it would stop counting. Cheapest change,
   but adds a second threshold to calibrate.
3. **Subtract an estimated pedestal measured from the negative controls.** Honest about the
   mechanism and uses data already produced every run, but the estimate is specimen- and
   background-dependent, so it is the least portable.

### Option 2 is now calibrated, 2026-10-05 — the neighbour's SHAPE separates the two cases

Option 2 ("require the neighbour to be credible before it can out-vote") carried the
objection that it "adds a second threshold to calibrate". That threshold has now been
measured rather than argued about, by `scripts/measure_neighbour_uniformity.py` over the
`runs/lod_v4` bundle BAMs. The coverage TSV carries `best_neighbour_reads` but no breadth for
the neighbour, so the quantity cannot be read off an existing artifact.

    uniformity = observed_breadth / breadth_expected_if_the_same_reads_were_uniform

A present organism's reads cover its genome, so uniformity approaches 1. Cross-mapped reads
pile onto the conserved fraction two genomes share, covering far less than their bulk
implies. It is a SHAPE statistic, near-independent of abundance, which is what lets one
threshold serve every spike level.

| bundle | declared neighbour | in the specimen? | neighbour uniformity |
|---|---|---|---|
| 632 *Y. pestis* | 630 *Y. enterocolitica* | **no** | **0.114 - 0.124** |
| 1392 *B. anthracis* | 1428 *B. thuringiensis* | yes | 0.639 - 0.973 |
| 1639 *L. monocytogenes* | 1642 *L. innocua* | yes | 0.935 - 0.992 |
| 487 *N. meningitidis* | 486 *N. lactamica* | yes | 0.970 - 0.979 |
| 28901 *S. enterica* | 562 *E. coli* | yes | 0.956 - 0.974 |

The gap runs from 0.124 to 0.639, so any threshold in between divides them, and the split
reproduces the scope this entry states independently: one defective bundle, four where the
gate does the job it claims. The 630 figure is constant across all four spike levels and
all three replicates, which is the pedestal signature - ~2,700 reads covering ~1% of a
4.6 Mb genome, about 8.5x more concentrated than uniform.

Why raw breadth cannot do this: breadth scales with abundance, so a trace true positive and
a fat pedestal can share one. At 100 rpm the target holds 2,336 reads and the absent
neighbour 2,731 - close enough that the raw ratio attributes the sample to the neighbour -
while their uniformities are 0.86 and 0.12.

⚠️ **The script's first verdict logic asked the wrong question** and is recorded because the
shape of the error matters: it tested "are the neighbours below the targets?", which assumes
the defective case. On a bundle whose declared neighbour IS present the ordering simply
inverts - the neighbour is credible and the unspiked target is the pedestal - so it printed
NOT SEPARATED for three bundles that are cleanly separated the other way round. The gate
under consideration asks only whether the NEIGHBOUR deserves to out-vote, so the neighbour's
own shape is the entire input; the comparison was a leftover from how the defect was first
seen.

**Not implemented.** This calibrates option 2; it does not choose it, and no gate has been
changed. The cliff measurement below is what confirms the mechanism the remedy would be
built on, and it is running as `d12fine`. Implementing a threshold before that lands would be
tuning against an unconfirmed model.

**SCOPE, measured 2026-10-03 and narrower than first stated.** The same read attribution
applied to the negative-control bundles shows the defect is NOT general. For five of the six
other taxa that stage 1 nominated, the declared near-neighbour is genuinely present in the
standing background - *B. thuringiensis* for *B. anthracis*, *L. innocua* for
*L. monocytogenes*, *E. coli* for *Salmonella* and for *Shigella*, *N. lactamica* for
*N. meningitidis* - so the denominator represents real competition and the gate does the job
it claims. That is why `attributed_to_neighbour` scored PPV 0.500 rather than 0: half those
calls were correct.

The defect needs BOTH conditions together:

  (a) the declared near-neighbour is ABSENT from the specimen, and
  (b) the neighbour reference recruits more background than the target reference does.

*Y. pestis* / *Y. enterocolitica* in this batch satisfies both. How often that pairing
arises across 89 taxa and real specimen backgrounds is unmeasured, and is the open question.

**CALIBRATION LIMIT, stated plainly.** Everything above comes from ONE simulated batch with
ONE background composition (2% *E. coli*, 1% *N. lactamica*, 0.7% each
*H. haemolyticus*/*L. innocua*, 0.5% *B. thuringiensis*) and reference pairs from one pack.
The 2,700/210 figures and the ~13x asymmetry are properties of those specific references and
that specific background. They should not be treated as constants.

**Would resolve it:** the finer rpm grid above (does the cliff land at ~180?), and -
for the calibration limit specifically - the TEN ZYMO ISOLATE RUNS, which are a better
instrument than the community arms for this question. Each isolate library contains ONE
organism, so recruitment onto any reference is unambiguous BY CONSTRUCTION, with no need for
the read-name trick that real community data cannot supply. Aligning each isolate against
each competitive bundle yields a per-reference recruitment coefficient in reads per million,
measured on REAL reads with real error profiles: a 10 x 89 matrix in place of the single
number this entry currently rests on. The even-vs-log standards then test whether recruitment
scales linearly with background abundance, and the *Salmonella* (neighbour present) versus
*Listeria* (neighbours absent) contrast tests condition (a) on real data.

What Zymo CANNOT settle: it contains no *Yersinia*, *B. anthracis*, *Francisella* or
*Burkholderia*, so those agents can be measured only as recruiting REFERENCES, never as
low-abundance targets; its ten organisms are a floor on background diversity, not a model of
a clinical specimen; and without the host-mixture arm it says nothing about the human-read
recruitment D13 records.

---

## D13 — Masking is detected at 200 bp / 95% identity but reads map at 150 bp; the two scales do not meet

**Finding:** `[MEASURED]` 2026-10-03, from the negative-control bundles of
`runs/lod_v2_kraken2` (no watchlist target spiked; reads carry their source accession).
Reads at MAPQ>=10 on each TARGET reference, by true origin:

| taxid | taxon | reads on target | dominant true origin |
|---|---|---|---|
| 5811 | *Toxoplasma gondii* | **8,310** | **human** NC_060925.1 / 060933.1 / 060943.1 |
| 487 | *Neisseria meningitidis* | 319 | *H. haemolyticus* (275), *E. coli* (42), *N. lactamica* (2) |
| 1392 | *Bacillus anthracis* | 66 | *L. innocua* (64) |
| 1639 | *Listeria monocytogenes* | 45 | *B. thuringiensis* (45) |
| 622 | *Shigella dysenteriae* | 61 | *E. coli* (60) |
| 28901 | *Salmonella enterica* | 22 | *H. haemolyticus* (20) |

Every cell is `unique_reads` from that sample's coverage TSV, i.e. MAPQ>=10 primary with
masked-region reads excluded - the quantity the pipeline actually gates on. Earlier versions
mixed units across the rows: 24,101 (any-MAPQ, with supplementary) for 5811, 317 for 487
(the sum of two of its three origins) and 23 for 28901 (a value in no replicate).

The *T. gondii* reference is a **human-read sink**: 8,318 host reads that survived depletion
map to it at MAPQ>=10, **26x** what the worst bacterial reference recruits (319 onto
*N. meningitidis*; 8,310/319 = 26.0). An earlier version of this table put 24,101 in that cell and claimed
"three orders of magnitude" - that figure was counted at ANY MAPQ including supplementary
alignments while every bacterial row beside it was MAPQ>=10 primary, so the comparison was a
unit artifact. All rows are now MAPQ>=10 primary. Its neighbour side holds **0** reads at MAPQ>=10 primary, so `neighbour_ratio` is
**`inf`** - which is what the pipeline's own TSVs record in 17 of 18 negative 5811 rows - and
the neighbour gate offers no protection whatsoever. An earlier version said "33 reads, so
730"; both halves were any-MAPQ figures left behind when the table cell beside them was
corrected to MAPQ>=10 primary. Breadth, occupancy and identity are
what actually hold the line - the specificity table shows 0.000 false positives - so this
is a redundancy that is doing nothing, not an active false positive.

**Why masking does not catch it.** `mask_pack_refs.py` runs at `min_hit_len=200` and
`min_identity=0.95`. A 150 bp read needs neither - and 150 bp is the length this batch
actually used (`simulate_benchmark.py:120`); an earlier version of this entry said 101 bp
throughout, which would have sent the remediation below to the wrong parameter: it maps at MAPQ>=10 on far shorter and
far less identical stretches than a 200 bp 95%-identity alignment. The result pack-wide:

* *T. gondii*: **1,982 bp** host-masked on a 65.7 Mb genome
* **2 of 359 genomes** have ANY host masking at all, 2,737 bp in total
* `crossgroup` is **0 bp for the two references this entry turns on** (630 and 632), and
  pack-wide only 20 of 359 genomes carry any at all, **81,196 bp** as the mask manifest
  totals it (19 genomes / 79,361 bp if the duplicate `620_Shigella` row is dropped -
  `README.md`'s convention; its BED is interval-identical to `622_Shigella_dysenteriae`,
  the same assembly under two watchlist rows, differing only in the header comment).

  The earlier figure of 81,565 bp was quoted from a manifest that disagreed with its own
  BEDs: `crossgroup_bp` read 1,835 for taxid 620 and 2,204 for 622 although their intervals
  are identical, because 620 had been masked standalone while 622 was carried forward from a
  shared-sequence run. Re-masking both pairs together (263/264 and 620/622) made the
  per-source columns agree and populated the `SHARED_SEQUENCE` flag that had been empty -
  the provenance marker `docs/STATUS.md` records as a known gap. `tests/test_pipeline.py`
  group `shared_assembly` now compares the per-source columns too; it previously compared
  only the four totals, all of which agreed, so it passed green on the inconsistency it is
  named for. It compares only across
  distant GROUPS by design, so bacteria-vs-bacteria recruitment is invisible to it - note
  *S. dysenteriae* (622), in the table above, does carry 2,204 bp of crossgroup from
  non-bacterial sources, so "0 bp for every genome checked" as this entry once read was
  wrong

So the masking step cannot see read-scale recruitment. This is a scale mismatch, not a
missing source.

**Consequence for D12 option 1.** Masking within-family recruitment was listed there as the
lever that would fix numerator and denominator together. At current parameters it would
not: the *E. coli* pedestal on *Y. enterocolitica* sits on a reference with 571 bp masked
and 0 bp crossgroup. Making masking catch it means read-scale parameters (~150 bp here,
and in general the batch's own read length, lower identity), which on closely related
genomes would mask large fractions of the reference and cost sensitivity in exactly the
silent way masking exists to prevent. That trade-off has to
be measured, not assumed.

**Decision: record, and do not change the masking parameters.** Both numbers that matter -
0.000 false-positive rate and 0.349% overall masked - are currently acceptable, and the
mechanism is now documented so a future false positive on a eukaryotic reference is
diagnosable rather than mysterious.

**Would resolve it:** per-reference recruitment measured on REAL single-organism reads, so
masking parameters can be chosen against a measured curve rather than a guess. The ten Zymo
isolate runs are exactly that instrument - see D12.

---

## D14 — A notifiable parasite is reported `supported` in every sample, from 100% human reads, and the benchmark does not count it ⚠️ *corrects an earlier "specificity is clean" claim*

**Finding:** `[MEASURED]` 2026-10-03, `runs/lod_v2_kraken2`. In a NEGATIVE CONTROL
containing no parasite at all, *Toxoplasma gondii* is reported at tier **`supported`**:

| metric | value | gate | verdict |
|---|---|---|---|
| unique_reads | 8,310 | floor | pass |
| breadth | **0.441** | >= 0.01 | pass |
| window_occupancy | **0.965** | >= 0.15 | pass |
| mean_identity | **0.9815** | >= 0.98 | pass |
| neighbour_ratio | **inf** | >= 1.5 | pass (no neighbour reads exist) |
| depth_cv | 5.44 | <= 2.50 | **FAIL** |

Every gate passes except `depth_cv`. That single threshold is the only thing between this
and `confirmed_in_silico` for a notifiable agent.

**The reads are 100% human.** Read attribution on the bundle BAM: of **8,318** records on the
*T. gondii* reference at MAPQ>=10 primary, **8,318 come from T2T human chromosomes**
(NC_0609*). Not a majority - all of them. (An earlier version of this paragraph said 8,380,
which is the `-F 0x4` count including supplementary alignments - the same unit error this
entry retracts further down. It was corrected in the retraction and left standing here.) 44% of a 65 Mb parasite genome is "covered" at 96.5% window
occupancy and 98.1% identity by host reads that survived depletion.

**It is systematic, not a *T. gondii* quirk.** Measured recruitment per million human reads
(`scripts/measure_recruitment.py`):

Measured at **150 bp** with **insert 350** and **MAPQ>=10**, matching `lod_v1`
(`scripts/simulate_benchmark.py:120,1202`; all 138 ART commands in that batch are `-l 150`),
over 10^6 simulated source reads:

| reference | records per 10^6 human reads | % of human reads |
|---|---|---|
| *Toxoplasma gondii* | **54,258** | **5.43%** |

The matrix now records `read_len`, `mapq_min` and `insert_mean` as columns, so a
coefficient can no longer be lifted off disk without its conditions.

⚠️ **This figure was published four times and wrong three times.** The history is kept
because the pattern is the lesson:

| | *T. gondii* | why it was wrong |
|---|---|---|
| first | 7.74% | SAM records divided by read PAIRS, labelled per-read - 2x high |
| second | 3.87% | unit fixed, but measured at 101 bp when the batch is 150 bp |
| third | 5.51% | a 50k-read spot check quoted while a 500k run was producing a better one |
| **now** | **5.43%** | 10^6 reads, 150 bp, insert 350, from the matrix file itself |

Every wrong version was hand-copied from a terminal into a document. The coefficient is
read-length dependent (1.4-2.0x between 101 and 150 bp, varying BY REFERENCE - it is not a
single factor), so `--read-len` is required, the artifact records it, and
`scripts/check_pedestal_model.py` now computes the comparison from primary sources and
REFUSES a matrix whose read length does not match the batch. Neither the ranking nor the
mechanism moved across the four versions, because those rest on direct record counts.

Eukaryotic references absorb host reads at percent-level rates. The *Zika* row deserves
separate attention: **0.26%** of human reads landing on a 10.8 kb RNA virus (2,555 records
per 10^6 at 150 bp) is not plausible homology and suggests human or vector sequence INSIDE
that reference - the Wright et al. mechanism, which would make it a contaminated reference
rather than a homology artefact. Flagged, not concluded.

**THE BENCHMARK IS BLIND TO IT, AND THIS CORRECTS AN EARLIER CLAIM.** All **18** `supported`
calls in the batch are *T. gondii*, one per sample, in LoD samples and negative controls
alike. The evaluator maps `supported` to POSITIVE - it is an assertion of presence. Yet
`specificity.tsv` contains **zero** *Toxoplasma* rows and the measured false-positive rate
was 0.000, because `negative_host_only`'s `truth.absent` lists eleven taxa and *T. gondii*
is not among them. An unlisted taxon is scored `unlisted`, not `false_positive`.

So "specificity is clean: 0.000 false-positive across all 6 negative controls for every
watchlist taxon", recorded in `docs/validation.md` on 2026-10-03, was wrong. It is 0.000
for the taxa the scenario happens to name. The most consequential false positive in the run
was invisible to the measurement - the project's recurring failure shape, this time in the
benchmark rather than the pipeline.

**What is actually holding the line.** Host-depletion efficiency is the control on this
class of false positive, and until 2026-10-03 nothing in the pipeline recorded it
(`ch_host_stats` was a dead channel; now written to `<outdir>/qc/host_depletion.tsv`). It
still GATES nothing.

⚠️ **Two earlier versions of this paragraph were wrong, the second in three ways.** It
first inverted the coefficient to claim "~152,000 surviving human reads", which exceeded
the sample's entire non-host remainder. The replacement then said "residual human reads are
**at most** 20,236, so the rate is ~41%, about 7.6x" - and that was wrong because:

1. **the inequality was inverted.** `residual = screened - nonhost_survivors`, and
   survivors are at most what was GENERATED, so the expression is a **floor** on residual
   host, not a ceiling;
2. the bound was **refuted by the same BAM it described** - 8,318 host-origin records land
   on this one reference, and the "at most" figure was for the whole sample;
3. the numerator **8,380 counted supplementary alignments** (`-F 0x4`) while the
   coefficient excludes them (`-F 0x904`), so the ratio mixed units.

`scripts/check_host_recruitment.py` now computes this, with the inequality direction in the
code and the filters forced to match:

    screened records           1,980,236
    residual host EXACT           24,708   host-origin primary records in the BAM
    records on ref 5811            8,318   MAPQ>=10, primary, all host-origin
    whole-genome coefficient        5.43%  (150 bp, insert 350, MAPQ>=10)
    in-pipeline rate EXACT          33.7%
    enrichment EXACT                 6.2x

A third version of this paragraph reported 35.1% / 6.5x as CEILINGS against a floor of
23,708, and its script asserted the exact count was unknowable because "the surviving reads
are not labelled". They are labelled - the numerator depends on that label - and the BAM
already open holds every screened read as a primary record, so one more pass gives the exact
denominator. The floor is still computed as a cross-check. The mechanism: survivors of host depletion are **not a random sample of host reads**. They are
precisely the reads that failed to map to CHM13 - repetitive, divergent, low-complexity -
which is the same class that cross-maps to a 65 Mb parasite genome.
`measure_recruitment.py` simulates from the whole host genome with no depletion step, so
its coefficient describes a population the pipeline never sees. **The coefficient must not
be inverted**, and "recruitment scales linearly with whatever survives" is wrong for the
host arm. It holds for *E. coli*, which passes depletion untouched - which is exactly why
that arm predicts within 1% while the host arm was out by 7x.

**The predictive model holds on both arms**, both computed by
`scripts/check_pedestal_model.py`. The neighbour arm covers the 12 samples with a 632
coverage TSV, and the script PRINTS its skips by reason: 6 negative controls that never
contained *Y. pestis*, and 3 `rpm1` samples that did but which stage 1 never nominated. An
earlier version of this paragraph asserted that 3/6 split while the script, testing the wrong
component field, printed all 9 under "stage 1 did not nominate it". The target arm is counted from the BAMs via `--attributed-bam`, because the
coverage TSV does not separate which of the target's reads came from *E. coli*; it replaced
a free-text flag that accepted bare integers with no filter or provenance.

    coefficient @150 bp, insert 350, MAPQ>=10, over 10^6 source reads
      E. coli -> ref 630 (neighbour)   3,355 rpm
      E. coli -> ref 632 (target)        256 rpm
    source reads   799,992  (component JSON reads_generated: 399,996 PAIRS x 2;
                              reads_requested is 400,000 - the plan, not the result)

    NEIGHBOUR ARM  predicted 2,684  vs observed 2,581-2,811   mean error  -0.9%  (n=12)
    TARGET ARM     predicted   205  vs observed   181-238      mean error  -2.9%  (n=12)

Both arms are counted by the script, under the coefficient's own `-F 0x904` filter, over the
SAME 12 samples. Two earlier versions of the target arm were wrong in different ways: -3.0%
from four hand-entered values measured at `-F 0x4` (including the supplementary alignments
the coefficient excludes), then **-0.6% from four of the twelve available BAMs with no record
of which four**. Fixing the filter while leaving the operator free to choose the samples moved
the laundering rather than removing it, so the script now DERIVES the BAM set from the samples
it just averaged - the two arms cannot silently differ in subject.

So `pedestal = recruitment_rpm x source_reads` is calibrated, within ~3% on both arms. The
target arm is reported deliberately: an earlier version quoted only the neighbour arm, which
is how a half-calibrated model gets called calibrated - and a later one quoted the target arm
from a subset, which is how it gets called better-calibrated than it is.

Four claims about this one number were published before this one, and three were wrong:
`276a3f1` "within 8%" (right by accident - a 2x unit error cancelling a 0.81x denominator
error), `d2d8fb9` "2.1x out" (wrong denominator and wrong read length), `f197d39`
"within 6%" (6.31% on its own two inputs, and both inputs were spot-check values superseded
by a larger run already on disk). The denominator must be the reads that ENTERED the
competitive alignment - the component JSON states it in PAIRS, so x2 - not a MAPQ-filtered
count of what survived.

**Decision: fix the BENCHMARK first, not the thresholds.** The measurement has to be able to
see the failure before any threshold is tuned against it, or the tuning is unfalsifiable:

1. `truth.absent` in every negative/specificity scenario must be the FULL watchlist minus
   what is deliberately present. Enumerating eleven taxa by hand guarantees the next gap.
   This is a one-line change in intent and a real change in coverage.
2. A `host_depletion_efficiency` figure belongs in the run metadata, because D14 shows it is
   load-bearing and it is currently unmeasured.
3. Only then consider `depth_cv`, which is carrying this alone, or host masking at read scale
   (D13), which would attack the cause rather than the symptom.

**Would resolve it:** re-running the batch after (1), which turns an invisible `unlisted`
into a counted false positive and gives the first honest specificity number for the pack.

### RESOLVED 2026-10-05 — the batch had no control that could contain the artefact

Step (1) landed and is confirmed working: the 2026-10-04 `divergent_v1` bake-off scores
*T. gondii* as `supported: 6 calls, 0 truth present, 6 truth absent, PPV 0.000`. It is
counted now, where it used to be `unlisted`. That made the real cause visible.

**The suppression machinery was never broken. It had nothing to compare against.**
`score_evidence.py` already has the right path - a taxon present in a negative control at a
comparable level is tiered `suppressed_control` - and it demonstrably works, suppressing
*E. coli* in all 18 samples of `lod_v3_threearm`. It never fired for *T. gondii* because
`build_ntc_plan` sets `host_fraction = 0.0`:

| | host reads | *T. gondii* finding |
|---|---|---|
| sample | 19,019,995 | 8,447 reads, 4,267 RPM, `supported` |
| `*_NTC` (role=control) | **0** | none possible |

An NTC is a *reagent* control: it models kit and index contamination. The dominant false
positive here is host-derived, so no quantity of NTC can ever contain it, `control_rpm`
stays 0 for taxid 5811, and the gate is unreachable. The pack had a no-template control and
no **matrix blank** - the negative extraction control that is the host matrix itself,
processed identically. In `divergent_v1` the gap is even starker: no control at all, and
the findings say so in their own reason string, `NO NEGATIVE CONTROL`.

**Demonstrated on data already on disk, before changing anything.** Re-running
`score_evidence.py` over `lod_v3_threearm`'s own coverage TSVs with the three
`negative_host_only_illumina_r[123]` samples relabelled `role=control` - they are already
host+background with no watchlist target, i.e. a matrix blank wearing the wrong label. The
offline baseline reproduces the pipeline exactly first (178 rows, **0** tier mismatches), so
the delta is attributable to the relabel alone:

| tier | before | after |
|---|---|---|
| `supported` | **18** | **0** |
| `suppressed_control` | 18 | 30 |
| `control_background` | 3 | 9 |
| `confirmed_in_silico` | **6** | **6** |
| `attributed_to_neighbour` | 74 | 74 |
| `screened_only` | 59 | 59 |

All 18 `supported` rows were *T. gondii*; 15 became `suppressed_control` and the 3 control
rows `control_background`. **Every *Y. pestis* call is untouched** - 6 `confirmed_in_silico`,
6 `attributed_to_neighbour` - so the fix costs no sensitivity on the one true positive the
batch contains. One matrix blank per BATCH is sufficient: `control_rpm` is built across the
whole run, not per scenario, which is why relabelling three samples suppressed the artefact
in twelve others from a different scenario.

**The fix:** `controls.include_matrix_blank`, implemented in
`scripts/simulate_benchmark.py:build_matrix_blank_plan` and enabled on `negative_host_only`.
It is the scenario's own host+background composition with no spike, `role=control`, at full
sample depth - the comparison is RPM against RPM and a shallow blank gives a noisy
denominator for exactly the trace artefacts it exists to catch. It does not replace the NTC
and does not make the negative scenario's own samples redundant: the NTC isolates reagent
contamination, the negative samples MEASURE the false-positive rate, and the blank is what
they are measured against. Asking for a blank on a scenario with no host is refused outright,
because that is an NTC under a different name - the guard was tested by constructing a
scenario that trips it, not assumed.

`tests/test_pipeline.py` group `matrix_blank` asserts BOTH halves: that a host-free NTC
leaves the artefact at `supported` (reproducing this entry's bug, so the test fails if the
control type is removed) and that a matrix blank moves it to `suppressed_control`.

**Confirmed end to end 2026-10-05 by `runs/lod_v4`.** The batch was regenerated with the
lod_v1 command (24 samples: 15 LoD + 3 negative + 3 NTC + 3 BLANK), and rep 1 of
`negative_host_only` regenerated at the original `--seed 7` is md5-identical after
decompression to lod_v1's, so the comparison is like for like. `supported` went 18 -> **0**,
`confirmed_in_silico` stayed at 6 with PPV **1.000**, the LoD stayed at **1000 RPM**, and
*T. gondii* now appears in `specificity.tsv` at `false_positive_rate` 0.000 with `worst_tier`
`suppressed_control` - counted and suppressed, where the table used to hold no Toxoplasma row
at all. The only suppressed taxa are *E. coli* (a third of the NTC) and *T. gondii*, so no
true positive was traded away. See docs/validation.md.

**Still open.**
`depth_cv` remains the only confirmation gate this artefact fails, and D13's read-scale
masking - which would attack the recruitment itself rather than suppress its product - is
still deferred pending the Zymo measurement. Suppression is the right control for a *known*
matrix; it does nothing for a specimen type whose matrix was never blanked.

---

## D15 — A divergent agent is invisible because it is classified CORRECTLY, and the watchlist registers its sibling as a suppressor ⚠️ *the amino-acid arm does not fix this, and was tested on the scenario built to favour it*

`divergent_bundibugyo_filovirus` and `divergent_clostridium_bont` exist to answer one question:
does the pipeline raise a lead for a divergent or novel agent, or does it go silent? Batch
`simulated/divergent_v1` (6 samples: Bundibugyo virus at 200 RPM x3, *C. butyricum* at 300
RPM x3) was run through `scripts/run_bakeoff.sh --variants "kraken2_minimizer plus_metabuli"`,
a true minimal pair over one samplesheet.

**The pipeline went silent, and the benchmark's own gate caught it.**
`divergent_bundibugyo_filovirus / Zaire ebolavirus / min_tier >= screened_only` -> **FAIL**,
`no finding` in all three replicates. No filovirus appears in `findings.tsv` for either
variant. Nor does any *Clostridium*.

**It is not a classification failure. Both classifiers are exactly right.**
Kraken2 and Metabuli each assign all **3,992 of 3,992** spiked reads to *Bundibugyo virus*
(taxid 565995) and roll them up identically through Filoviridae -> Orthoebolavirus ->
*Orthoebolavirus bundibugyoense*. Not a single read differs between the two arms. Zaire and
Sudan ebolavirus receive **zero** reads from either.

**The mechanism is the watchlist's shape, and it is the opposite of a gap in sensitivity.**

  - The resolved watchlist holds 89 entries, **71 `species` and 18 `no rank` - nothing above
    species**. There is no Filoviridae or Orthoebolavirus entry for a divergent filovirus to
    trigger.
  - Zaire ebolavirus (186538) lists `565995` - Bundibugyo - in `near_neighbour_taxids`. The
    pipeline KNOWS the relationship. But near neighbours are the competitors bundled into the
    stage-2 index, the organisms a target must out-recruit to earn `neighbour_ratio`. They are
    **suppressors, not triggers**.

So reads from a divergent sibling are classified correctly to a non-watchlist taxon, and the
only role that taxon can ever play is to argue AGAINST a watchlist call. The better the
classifier, the more completely the agent disappears. A tool that was WRONG - that smeared
those reads onto Zaire - would have raised the lead the scenario expects.

**That also makes the scenario's expectation unsatisfiable as written.** `min_tier:
screened_only` on *Zaire ebolavirus* presumes detection by cross-mapping onto a sibling, which
only occurs when the classifier misassigns. The expectation encodes the wrong mechanism, so it
cannot be met by a correct pipeline. It is kept failing rather than relaxed, because the
underlying gap is real: the thing it is pointing at should be fixed, not the assertion.

**The amino-acid arm does not address this, and this was its best case.**
`findings.tsv` is identical between the two variants: same 48 rows, `tier` differing in **0**
of them, only `n_classifiers`/`classifiers`/`screen_reads`/`reasons`/
`mapped_fraction_of_assigned` moving. Evaluator totals are identical (0 TP, 0 FN, 6 FP, 0
false confirmations). Cost: 153.2 -> 173.7 CPU-minutes for 6 samples. Combined with the
identity-matched result in `docs/benchmarking.md`, Metabuli has now changed **zero decisions
on every scenario it has been run against**, including the `divergent_*` pair the variant
table nominated as where its value "should show up". Its remaining untested case is an agent
absent from the index, which no scenario in `human_core.yaml` constructs, and which cannot be
built from a database (see the DIVERGENT AGENT PROXY comment).

**Options, none chosen here.**

  1. **Rank-aware watchlist entries.** Allow a genus/family row (Orthoebolavirus, Filoviridae)
     at lead tiers only, so a sibling assignment raises `screened_only` without ever reaching
     a confirmation tier. Costs specificity on commensal-rich families; the negative-control
     scenario sets the price.
  2. **Promote near neighbours to triggers at lead tier.** `near_neighbour_taxids` already
     encodes the relationship; a read mass on a near neighbour of a high-consequence taxon
     could raise a flagged lead naming the neighbour. Keeps the species call honest, but
     doubles as a second use for a field that currently means "suppressor", which invites
     exactly the kind of unit confusion D12 and D13 came from.
  3. **Lean on the marker/HMM tier instead**, which is what the scenario text actually
     anticipates ("Expect at least a flagged lead via the HMM tier"). `MARKER_SEARCH_HMM` ran
     on these samples and produced no filovirus lead either; why is unresolved and is the
     cheapest thing to investigate first, since it needs no watchlist change.

VERIFY: whether `MARKER_SEARCH_HMM` has filovirus profiles at all, and if so why 3,992
on-target reads produced no marker evidence. That is the next question, and it is independent
of the classifier arms.

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
