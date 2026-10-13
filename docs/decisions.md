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

### ⚠️ REMEDY DEFAULTED OFF 2026-10-08 — real data refuted its threshold, see D18

Everything below is the simulated result and it still stands AS A SIMULATED RESULT: the
cliff really does disappear on `d12fine`. What does not transfer is the threshold. On the
Zymo libraries the gate took false positives from 4 to 25 on the isolates and 0 to 4 on the
community, and did not rescue the true positives it exists for, because their neighbours
scored ABOVE 0.50. `min_neighbour_uniformity` now defaults to 0. The undeclared read floor
this entry documents is REAL and still unfixed; what is retracted is that this particular
threshold fixes it.

### REMEDY IMPLEMENTED 2026-10-07 — option 2, and the cliff is gone

Option 2 ("require the neighbour to be credible before it can out-vote") is implemented. A
neighbour may only win the attribution if its reads COVER its genome rather than piling on
sequence it shares with the target:

    uniformity = covered_bp / min(reference_bp, summed_reference_span)

`coverage_stats.py` emits it per finding as `best_neighbour_uniformity`; `score_evidence.py`
gates on `--min-neighbour-uniformity` (default 0.50, `params.min_neighbour_uniformity`).
Absent evidence counts as CREDIBLE — no neighbour, or samtools coverage unavailable, leaves
the gate exactly as it was — so the change can only ever withhold an attribution on positive
evidence of a pedestal, never manufacture one.

**The cliff this entry predicted and then measured at 176.6 rpm is gone.** Re-scoring
`runs/d12fine`, where the declared neighbour is absent from every sample and scores
uniformity 0.146–0.159 at every level:

| level | before | after |
|---|---|---|
| 120 / 150 / 170 | 3/3 `attributed_to_neighbour` | **3/3 `supported`** |
| 190 / 210 / 250 | 3/3 `supported` | 3/3 `supported` |

*Y. pestis* is now judged on its own evidence at 120 rpm, the lowest level tested, instead of
being attributed to an organism that is not in the specimen.

**The cost, measured on `runs/lod_v4` rather than assumed.** This is a PERMISSIVE change, so
both columns belong in the record:

| | before | after |
|---|---|---|
| `attributed_to_neighbour` | 84 | 77 |
| `supported` | 0 | 7 |
| `confirmed_in_silico` | 6 | 6 |
| actionable TP / FP | 6 / 0 | **12 / 1** |

Six new true positives — *Y. pestis* at 10 and 100 rpm, genuinely spiked and previously
attributed away — against one new false positive, *P. vivax*. That is the trade, and it is
the trade this entry argued for: the attribution was wrong in exactly the case where the
neighbour is not present, and removing it recovers detections at the cost of one call that
the confirmation gates then have to carry.

⚠️ **TWO THINGS I GOT WRONG, BOTH CAUGHT BY THE VALIDATION RATHER THAN BY REVIEW.**

1. **The threshold's justification was weaker than I first stated.** It was calibrated on 5
   bundles, which gave 0.1532 (absent) against 0.9720–0.9913 (present) and made 0.50 look
   like the middle of a wide empty gap. Across all 182 measured rows in `lod_v4` the nearest
   values straddling the threshold are **0.4321 and 0.5010** — a 0.07 gap. The distribution
   is still strongly bimodal (19 rows below 0.3, 162 above 0.5), and both straddling rows are
   controls, so no actionable call in this batch turns on where exactly 0.50 falls. But it is
   not the clean separation the first calibration implied, and a batch with a different
   background could place a real decision in that gap.

2. **The first implementation did not do what its own comment said.** The not-credible case
   sat in its own `elif`, which appended the reason and then skipped the gate block, leaving
   `tier` at its initial `screened_only`. The comment claimed it "falls through to the
   confirmation gates"; it did not. The symptom was that the d12fine cliff did NOT move —
   120/150/170 simply became `screened_only` instead of `attributed_to_neighbour`, so a 170
   rpm sample scored below a 190 rpm one not because its evidence was weaker but because it
   was never evaluated. Only re-running the cliff measurement exposed it; the unit tests
   passed throughout.

**Still open.** The calibration rests on one simulated batch with one background composition,
the same limit this entry records for everything else in it. The Zymo isolate runs remain the
instrument that would replace it with a per-reference measurement on real reads.

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

⚠️ **CORRECTED 2026-10-06: the *T. gondii* row's coverage geometry above is inflated
~1000x.** `samtools depth -a` omitted every zero-coverage contig, so breadth and occupancy
were computed over only the contigs that recruited reads - 57,470 of 65,668,120 positions on
that 2,264-contig reference. The READ COUNTS in the table (8,310 on target, the 26x ratio
against *N. meningitidis*) are from `collect_alignments`, not from depth, and are unaffected;
the breadth and occupancy figures quoted anywhere from this entry are not. Fixed in
`bin/coverage_stats.py` with `-aa` restricted to the target's contigs; `runs/lod_v4` re-scored
2026-10-06 gives breadth 0.000336 and occupancy 0.001256 for the same sample. The conclusion
of this entry - that masking operates at a scale reads do not reach - is unchanged, and the
*T. gondii* reference is a worse human-read sink than the entry implies, not a better one.

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

### CORRECTED 2026-10-06 — the breadth figures above are inflated ~1000x, and the diagnosis was wrong twice

**The numbers in the table at the top of this entry are wrong.** `bin/coverage_stats.py` ran
`samtools depth -a`, which emits positions only for references that received a counted read,
so breadth and occupancy were computed over only the contigs that recruited something. On
this exact bundle - 2,264 target contigs, 65,668,120 bp - `-a` emitted 57,470 positions,
0.087% of the reference. Corrected with `-aa` restricted to the target's contigs
(`runs/lod_v4`, re-scored 2026-10-06):

| metric | as recorded above | corrected |
|---|---|---|
| breadth | 0.441 | **0.000336** |
| window_occupancy | 0.965 | **0.001256** |
| mean_depth | - | **0.0165x** |
| depth_cv | 5.44 | **198.79** |

So "44% of a 65 Mb parasite genome is covered at 96.5% window occupancy" was never true. The
truth is starker and simpler: 8,318 host reads pile onto 0.03% of the genome. The mechanism
this entry describes is unchanged and the read-attribution evidence stands; only the coverage
geometry was mismeasured.

**The sentence "Every gate passes except `depth_cv`" was wrong, and so was my replacement for
it.** With corrected depth, breadth and occupancy fail too - the reasons string reads
`failed: breadth 0.0003 >= 0.0050; occupancy 0.0013 >= 0.1500; depth_cv 198.79 <= 2.50`. I
predicted that this would demote the finding to `screened_only` and that the gates would
therefore have caught the false positive on their own. **It does not.** Re-scoring `lod_v4`
with the matrix blanks demoted back to samples, with the corrected metrics, still produces
*T. gondii* at `supported` in all 21 rows.

**That is the finding this correction exposes, and it is worse than the one it replaces.**
`supported` is not a judgement, it is the FALLBACK: `score_evidence.py` assigns it to
anything that clears `unique_reads` and `neighbour_ratio`, however badly it then fails the
confirmation gates. A finding that fails EVERY gate - breadth by 15x, occupancy by 119x,
depth_cv by 80x - is reported at the same tier as one that fails a single gate narrowly, and
the evaluator maps `supported` to POSITIVE. The tier cannot distinguish "nearly confirmed"
from "comprehensively refuted".

So the matrix blank is not belt-and-braces over a gate that would have held; **it is the only
thing standing between this artefact and a positive call**, before and after the depth fix.
That raises the value of the control recorded below and lowers the confidence this entry
should place in the gate structure.

✅ ANSWERED 2026-10-08 ON REAL DATA — see D18. In the Zymo single-organism libraries, where
every watchlist call but the library's own organism is false by construction, true positives
fail a median of 1 confirmation gate and false positives 5 (median breadth 0.8661 against
0.0015). Requiring breadth AND occupancy to PASS keeps 4 of 4 true positives and cuts false
positives from 45 to 4; `--actionable-requires` implements it, default `breadth,occupancy`.
Neither simulated benchmark moved. The original question, kept because it is what was asked:

VERIFY: whether a finding that fails N gates should fall to `screened_only` rather than
`supported`, and what N is. That is a change to the tier ladder, not a threshold, so it wants
its own decision and its own measurement of what it would reclassify across the 89 taxa.

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

## D20 — Three adversarial passes over everything, and what calibrating on two datasets actually cost

2026-10-13. Three independent reviews — threshold provenance, the open backlog, and an
attack on D14/D18/D19 — run against the code and the data on disk, not against the
comments. Every figure below that is marked ✓ I re-derived myself; the rest is the
reviewer's measurement, attributed as such.

### The class of defect the validation set cannot show

Not "the thresholds are loosely calibrated". Something sharper: **a gate that is
unreachable on realistic input is invisible on everything this project has ever tested.**

| gate | unreachable when | evidence |
|---|---|---|
| `min_alignment_length: 100` | any 2×100 library | ✓ **81 of 82** real Zymo rows fail it |
| `max_depth_cv: 2.5` | mean depth < 0.26× | no row below 0.258× has passed, in 709 rows |
| `min_read_identity` illumina 0.98 | strain-divergent real isolate | *L. monocytogenes*, 12.5M reads, markers PRESENT, held at 0.9643 |
| `min_occupancy: 0.15` | large genomes | *derived, not measured*: ~20,000 reads for *T. gondii* (✓ 66,891 windows, so 0.15 needs 10,034 occupied; observed 82 on 8,261 reads), and *S. mansoni*, which has no coverage row on disk |

✓ **Every `confirmed_in_silico` call in the project's history is *Y. pestis* or *E. coli*
on simulated 150 bp reads.** ✓ 13 runs carry a findings.tsv and 10 contain at least
one confirmation; not one is on real data, and the reason is the first row of that table
rather than anything biological.

✓ The unit suite cannot see any of them: the only `mean_alignment_length` fixture is set
to exactly `100.0`, the one value that passes. The standing lesson "a threshold must be
reachable, not merely strict" (D10, D11) is recurring for the fourth time because **the
tests check that gates reject bad data and never that they accept good data.** The
generalisable fix is a parameterised reachability test — for each gate, assert some
plausible real input passes, swept over read length 100/150/250 and both platforms.

Related, and the same shape: ✓ all 89 watchlist rows carry explicit `min_unique_reads`
and `min_breadth`, so the pack-level defaults are unreachable, and ✓ the identity
resolver always returns from the platform block, so `min_read_identity` and
`min_read_identity_by_group` are dead — while `nextflow.config` documents all of them as
live overrides. An override that prints a confirmation line to stderr and changes nothing
is worse than no override.

### The uniformity metric: a retraction of my own correction, and what is actually wrong

I reported that `coverage_stats.py` computed an unnormalised quantity while
`scripts/measure_neighbour_uniformity.py` computed a normalised one, and that D12's 0.50,
D18's retraction and D19's 0.30 were decided across that gap. **That was wrong.** `span`
is the summed ALIGNED reference span — the normalising term itself — so the two
expressions are the same formula: ✓ over all 467 neighbour rows, max |old − new| =
**0.0073**, mean 0.00015. The docstring I called false was correct.

What is real, and was the actual finding: **normalising never made this
abundance-independent, and the docstring claimed it did.** ✓ One ABSENT neighbour
(*L. innocua*), its own read count against its uniformity:

    14 reads 0.968 | 1,083 reads 0.852 | 11,593 reads 0.418 | 190,739 reads 0.265
    ...then 1,880,316 reads 0.364 and 7,639,529 reads 0.391 - IT RISES AGAIN AT THE TOP.

The decline is NOT monotone, and this entry first presented it as though it were by
stopping at the fourth point. Two further caveats on the same series: the first two points
sit below the 0.05x depth floor this entry itself calls uninformative, and the points come
from different specimens, so abundance is confounded with composition. The same organism
at 276,000 reads scores 0.9915 in `lod_v4` - higher than at 1,083 reads - because there it
is genuinely PRESENT. The mechanism claim survives all of that; the series as first
presented was selective.

Below ~1× depth the statistic is covered-per-aligned-base, which is ~1 for any taxon too
sparse for its reads to overlap; above ~1× it saturates to breadth, which for a sibling
plateaus at the shared fraction. It is a **depth-gated pile-up detector** described as a
level-independent shape statistic.

### The recalibration, and why it still does not license a threshold

All 512 competitive BAMs re-measured, truth labelled per specimen in three classes —
because the sibling case is the one that matters and lumping it into "absent" is what
made every previous calibration look clean. ✓ Counting DISTINCT CONDITIONS, not rows,
since 288 of the 292 present rows are one simulated background replicated:

| class | conditions | range (at ≥0.05× depth) |
|---|---|---|
| PRESENT | 9 | 0.524 – 0.991 |
| ABSENT_SIBLING | 17 | 0.022 – 0.861 |
| ABSENT_DISTANT | 26 | 0.000 – 0.604 |

⚠️ **Two caveats the first version of this entry left implicit.** The table is
DEPTH-FILTERED at ≥0.05×; unfiltered, ABSENT_DISTANT runs 0.000–1.000, so its 0.604
ceiling exists only because the floor presented below as *conclusion 1* was already
applied. And the three counts are not one unit: PRESENT collapses 288 replicate rows of
one simulated background to 5 conditions and adds 4 individual real rows, while the other
two count distinct (run, sample-family, neighbour); under the rule giving 17 and 26,
PRESENT is 19. The PRESENT arm is thinner than "9" suggests — **5 conditions are one
synthetic background and all 4 real-data conditions are the same organism, *E. coli*.**

1. ✓ Below ~0.05× neighbour depth the number is noise: unfiltered, ABSENT_DISTANT reaches
   p90 **0.932** from sparse non-overlapping reads alone.
2. ✓ Above that floor the arms appear to cross — present *E. coli* in `zymo_log` is
   **0.524**, absent *Citrobacter* in `zymo_even` ONT is **0.604** — but THAT PAIR IS
   WEAKER EVIDENCE THAN IT LOOKS and this entry first overstated it. The two points come
   from different platforms, and *Citrobacter* is probably mislabelled: `truth.py`'s
   SIBLING_OF list has no genus-level Enterobacteriaceae entry, so *Citrobacter*,
   *Klebsiella* and *Enterobacter* fall to ABSENT_DISTANT however close they are. At 72.7×
   depth inside the *E. coli* bundle of a specimen that contains *E. coli*, it is a
   sibling by behaviour; relabelled, the cross disappears. What does NOT depend on it: a
   0.63 threshold still scores **99.7% on rows** while calling present *E. coli* a
   pedestal, and conclusion 3 is sufficient on its own to keep the split unshipped.
3. ✓ ABSENT_SIBLING spans the whole range and always will. A sibling reference IS covered
   by its relative's reads. That is the honest answer and the reason the tier exists.

**`params.contest_attribution` is OFF.** The tier is kept; the pedestal/sibling split is
not shipped. Recommended resolution: drop the split rather than re-threshold it — a
contested row becomes `indistinguishable_complex`, which needs no number — but that is
coupled to whether the tier should rank above `supported`, which is unresolved.

Also fixed: a missing measurement now leaves the decision unchanged instead of
manufacturing a tier. ✓ `runs/lod_v1` and `runs/lod_v3_threearm` have no
`best_neighbour_uniformity` column, so byte-identical rows scored differently across
runs, and the lower-ranking tier would have failed the very expectation D19 cited.

Published, as the measurements D12 and D13 both name and neither produced:
`benchmarks/measurements/cross_recruitment_matrix.tsv` (1,773 per-taxid rows on real
reads), `neighbour_uniformity_recal.tsv` (467 rows, three-class truth) and
`neighbour_truth_labels.py`, which is the RULE that produced the `label` column. The
first version published the labels without the rule, which makes a derived column
unauditable - and the rule is hand-built, single-sourced and has at least one known
limitation (no genus-level Enterobacteriaceae sibling entry, see conclusion 2), so it is
the part most in need of review.

### What the alignment-length gate was actually worth

Sweeping `min_alignment_length` from 100 down to 0 across five batches and 512 rows
changes **exactly one call, and it is a true positive**. No row anywhere on disk with
`mean_alignment_length < 100` passes the other gates. The gate's measured record in this
project is **zero false positives prevented, one true positive suppressed** - a stronger
argument for fixing it than unreachability, and the one this entry should have led with.
It also means no value in [0, 98.9) is distinguishable from any other on current data, so
the Illumina setting is defended by read-length reasoning, not by measurement.

**AND THE TWO CHANGES TOGETHER NET TO ZERO NEW CONFIRMATIONS ON REAL DATA.** The
reachability fix confirms *Salmonella enterica*; D17a immediately caps it at `supported`.
`zymo_rl100` goes from `conf=0, sup=5` to `conf=0, sup=5`. "Produces the first real-data
confirmation in project history" is true of the first change in isolation and FALSE of
what ships, and in a project with this record of writing down justifications for
behaviour the code does not have, that distinction has to be on the page. The only
evidence the reachability fix works is the new synthetic fixture. What changed is the
reason attached to the row: from "failed: mean_alignment_length 99 >= 100", which was
false, to "all confirmation gates passed ... CAPPED AT SUPPORTED: no marker rule block",
which is true and actionable.

### Backlog, verified against current code

**Before any real use.** ✓ `min_alignment_length` (above). ✓ `--threshold
min_read_identity` silently dead. **D17a**, which is dormant only because nothing can
currently confirm: ✓ re-scoring at a reachable alignment floor produces the first
real-data confirmation in project history, and it is *Salmonella enterica*, `notifiable`,
`marker_verdict` EMPTY, with no block in `marker_rules.yaml` — so the two must land
together. **D16a** is no longer a prediction: injecting one 1,000 RPM control row into
`zymo_rl100` flips **2** genuine `supported` detections to `suppressed_control` (both the
*L. monocytogenes* trace contamination D18 established as real) plus 3 `screened_only`
rows — a true-positive → false-negative flip. `d12fine` carries 9 *Y. pestis* `supported`
rows at 2,149–2,925 RPM, every one of which a 1,000 RPM control would suppress at
`control_ratio 5.0`. The mechanism is sharpest when the operator does the right thing:
adding a control simultaneously removes the "NO NEGATIVE CONTROL IN RUN" caveat and
suppresses the finding.

**Before release.** **D15's cause is structural and was misdiagnosed**: the filovirus
markers DO fire (`ebola_NP`, `ebola_GP`, `ebolaS_NP`, `ebolaS_GP` all PRESENT, both
`confirm_requires` satisfied) and are then discarded, because `score_evidence.py` builds
findings only from coverage rows and the `--missing-refs` list — **marker evidence can
gate a finding but can never raise one**. The entry's "why is unresolved" is answered and
its VERIFY ("does HMM have filovirus profiles at all") is resolved YES. The fix is a
capped `marker_lead` tier, and it must land with **D17b**, since the same data shows
`marburg_L` PRESENT on HMM bits alone across a genus boundary — the HMM branch never
checks the catalogue's `min_coverage`/`min_identity`. Then **D16b** (units and threshold
together), and **declaring** D12's residual floor, which D19 moved rather than removed:
it is now ≈55–58 rpm and set by `min_occupancy`, i.e. by genome size rather than by
background composition. Predictable is better than unpredictable; undeclared is still
undeclared, and `min_occupancy` is a hardcoded argparse default with no per-taxon
override and no representation in the pack.

**Accept and document.** D13 (do not touch masking; the recruitment matrix it wanted is
now published), the Illumina identity floor as a strain-divergence problem rather than an
error-rate one, D17c.

**Explicitly rejected.** Replacing `min_occupancy` with an occupancy *ratio*: measured, it
buys 3 simulated true positives for 3 real-data false positives.

### What this says about the evidence base

✓ Across 709 stage-2 rows in six runs: 623 bacteria, 69 parasite, 17 fungus, **0 virus**.
30 of 89 watchlist entries are viruses — Variola, the ebolaviruses, Nipah, MERS — and not
one has ever produced a coverage row, so every virus-specific threshold is unexercised
rather than weakly calibrated. ONT has 11 rows, all from the batch that never completed.

And the thresholds have never been tested where they will be used: 65% of all rows sit in
depth regions where no placement matters (47% is cross-mapping noise, 18% passes
everything) and **2% fall in the band where a gate actually decides**. A clinical specimen
lives in the 0.01–1× band, which is ✓ **34.6%** of this corpus — corrected from
8.5%, which is the 0.05–1× band. The corpus is better populated there than first stated;
what remains true is that the deciding band is 2%.

---

## D19 — The attribution gate never looked at the target, and a 100%-abundance notifiable pathogen produced no finding at all

D18 item 2 recorded that near-identical pairs are silently attributed. Measuring what the
run actually *reported* turned that from an accuracy complaint into a structural one. In
the Zymo *C. neoformans* isolate library - one organism, 5,870,155 reads on its own
reference, breadth 0.642 - the complete output is:

| tier | taxon | reads | breadth |
|---|---|---|---|
| `attributed_to_neighbour` | *Cryptococcus neoformans* | 5,870,155 | 0.642 |
| `attributed_to_neighbour` | *Cryptococcus gattii* | 8,530 | 0.005 |
| `screened_only` | *Salmonella enterica* | 37 | 0.001 |

Nothing is reported. *C. deneoformans*, which took the vote at ratio 0.7072, is a
near-neighbour reference and **not a watchlist taxon**, so it generates no finding of its
own. The reads leave the pipeline unaccounted for. The *E. coli* library behaves the same
way via *Shigella* at 0.9032. Both organisms are notifiable.

**The cause is that the attribution branch never read the target's own evidence.** It
compared two read counts and a ratio. A target covering two thirds of its reference lost
to a sibling on a 10% count difference, and the loss was recorded in prose an analyst
could not even decode, because the neighbour was named by taxid.

### Two changes, and the order matters

**1. Name the neighbour.** The finding now reads `reads favour Shigella [taxid 620]
(near-neighbour reference)`, from the pack's existing refs manifest (`--refs-manifest`;
absent degrades to the bare number and changes no tier). `best_neighbour_name` travels as
a column. 45/45 attributed rows in `zymo_rl100` carry a name. This makes the call
auditable; it does not make it right.

**2. `--contest-attribution` (`params.contest_attribution`, default on).** Where the
TARGET's own breadth and occupancy pass - D14's gates, chosen because D14 measured them as
the discriminator - a one-sided attribution is no longer allowed to stand silently. The
neighbour's own coverage uniformity then decides which of two different situations it is:

| neighbour uniformity | meaning | outcome |
|---|---|---|
| ≤ 0.30 | reads pile on shared sequence; not in the specimen | attribution rejected, target judged on its own evidence |
| > 0.30, or unmeasured | both references genuinely covered | **`indistinguishable_complex`** |

`indistinguishable_complex` is a new tier below `supported` and above
`attributed_to_neighbour`: *one of this pair is present and this evidence does not say
which*. It is the honest statement, and D18 argued for exactly it.

### Why this threshold is not the one D18 retracted

It is the same metric and the same 0.30-ish region, so the difference has to be stated
precisely. `min_neighbour_uniformity` asked "may this neighbour out-vote?" of **every**
row, and on real reads let 25 cross-mapping artefacts through. This asks a narrower
question of a far smaller population - where the target's own coverage geometry passes,
which is 6 of 45 rows in `zymo_rl100` and **0 of 84** in `lod_v4`. In that population the
separation is not the 0.07 gap D18 refused to trust:

    simulated, neighbour absent   0.1495 - 0.1532   (12 rows)
    real, neighbour a sibling     0.4107 - 0.8607   (10 rows)

All 22 measured rows sit outside a 0.26-wide empty gap, and 0.30 is independently the
`--pedestal-max` already chosen in `scripts/measure_neighbour_uniformity.py`. The guard
that does the real work is the precondition, not the threshold.

### Measured, off against on, over four batches

The off column reproduces D18's recorded scoring **identically** in all four, which is
what makes the comparison worth reading.

| batch | `attributed_to_neighbour` | `indistinguishable_complex` | `supported` | false positives |
|---|---|---|---|---|
| zymo_rl100 | 45 → 39 | 0 → 6 | 5 → 5 | **2 → 2** |
| zymo_rl150 | 7 → 3 | 0 → 4 | 1 → 1 | **0 → 0** |
| lod_v4 | 84 → 81 | 0 → 0 | 0 → 3 | **0 → 0** |
| d12fine | 84 → 75 | 0 → 0 | 9 → 18 | **0 → 0** |

`screened_only`, `suppressed_control`, `control_background` and `confirmed_in_silico` are
unchanged in every batch; no row moved anywhere except out of `attributed_to_neighbour`.

**D12's defect is closed, and the evaluator says so rather than me.** `lod_v4` off FAILS
its own declared expectation - `lod_gram_negative/Yersinia pestis min_tier_at_rpm[100]:
expected >= supported, observed attributed_to_neighbour` - and passes with the flag on.
Twelve genuine *Y. pestis* spikes across the two simulated batches were being attributed to
*Yersinia* at uniformity 0.15, an organism absent from the specimen. That is the defect
D12 described; the remedy it proposed was too broad, and the missing piece was D14's gates
as a precondition.

### `indistinguishable_complex` COUNTS AS A DETECTION — decided, with the cost stated

Asserting that one of a named pair is present is an assertion of presence, and for a
notifiable organism a scorecard that calls it a miss is wrong about what the tool did. It
is in `POSITIVE_TIERS`; it is deliberately NOT in `CONFIRMED_TIERS`, because it resolves
no organism and must never read as a confirmation.

The resulting scorecard, with the before column reproducing D18's recorded figures
exactly:

| batch | before TP/FP | after TP/FP |
|---|---|---|
| zymo_rl100 | 3 / 2 | **7 / 4** |
| zymo_rl150 | 1 / 0 | **4 / 1** |
| lod_v4 | 6 / 0 | **9 / 0** |
| d12fine | 9 / 0 | **18 / 0** |

**Every one of the three new false positives is *Shigella dysenteriae*** in a specimen
whose *E. coli* is real, and that is not a scoring artefact to tune away. *Shigella* is a
pathotype within *E. coli*; the watchlist's own note says species-level separation from a
metagenome is unreliable. A method that cannot separate the pair should be charged for
both sides of it, and the number is the honest size of the problem `ipaH` plus lineage
markers (D17) would solve. The other two in `zymo_rl100` are the *L. monocytogenes* rows,
unchanged by this and already established as **real trace contamination** rather than
errors - confirmed by Kraken2 at 0.28% and 0.06% with breadth matching coverage theory, so
the truth model, not the pipeline, is what is wrong about them. Counting only calls the
tool actually got wrong, `zymo_rl100` reads 7 TP / 2 FP and `zymo_rl150` 4 TP / 1 FP, and
every remaining false positive in both is the *Shigella*/*E. coli* pair.

The remaining expectation failures on the Zymo batches are `zymo_even`'s **ONT** rows,
which have no findings because that batch was never completed - not a scoring outcome.

### What is NOT decided

**The pairs themselves are still unseparated.** This reports the limit instead of hiding
it. `ipaH` plus lineage markers for *Shigella*, and a *C. gattii*/*deneoformans*
discriminator, are what would actually resolve them - which is D17's territory, not this
entry's.

---

## D18 — The first real sequencing data, and it refuted one change and justified another

`real/mock_v1`, fetched 2026-10-07: the ZymoBIOMICS D6300 standard as EVEN and LOG
communities on both Illumina and ONT, plus **all ten member organisms as individual
Illumina libraries**. 26 files, 54 GB, every one md5-verified against ENA, no run refused
for zero reads, measured read lengths matching the declared `bracken_readlen`.

Run as three batches because Bracken takes ONE read length per run and
`modules/local/detect.nf` passes a single `-r ${params.bracken_readlen}`: `zymo_rl100`
(ten isolates + LOG community, ~2x101), `zymo_rl150` (EVEN community, ~2x151), `zymo_ont`
(both communities, GridION). A mixed run would not error - it would produce plausible
abundances that are simply wrong for the mismatched group.

**WHY THE ISOLATES ARE THE INSTRUMENT.** Each library contains exactly ONE organism, so
every watchlist call but that organism is false BY CONSTRUCTION. There is no simulation
assumption behind the truth, and no read-name trick is needed - which is what D12 means
when it names these runs as the answer to its calibration limit.

### What they refuted

**D12's neighbour-credibility remedy, measured the day before on simulated data.** On
`lod_v4` it read +6 true positives for +1 false positive. On real reads, in both
directions:

| | FP with gate | FP without |
|---|---|---|
| zymo_rl100 | **25** | 4 |
| zymo_rl150 | **4** (incl. *B. anthracis*) | 0 |

and it did not rescue the true positives it exists for: *Salmonella* at 2,067,999 reads
and *E. coli* at 541,547 were still attributed away, their neighbours scoring 0.7473 and
0.6172 - above the 0.50 threshold, hence "credible".

The mechanism is the interaction with the tier ladder, not the metric. Withholding a bad
attribution is correct; those neighbours really are absent. But withholding is NOT evidence
the target is present, and `supported` was the tier a finding reached by FAILING. Every
false positive gained in `zymo_rl150` had neighbour uniformity below 0.50. `min_neighbour_uniformity`
now defaults to 0; the measurement and the `best_neighbour_uniformity` column are kept,
because the separation is real on the bundles it was measured on - what is not established
is that the threshold transfers to a real background.

### What they justified — D14's VERIFY, answered

D14 asked whether a finding failing N confirmation gates should fall to `screened_only`
rather than `supported`. In the isolates, true positives fail a median of **1** gate and
false positives **5**; median breadth **0.8661** against **0.0015**, median depth_cv **0.7**
against **78.9**. The gates were never wrong - the ladder discarded their answer.

| rule | TP kept | FP kept |
|---|---|---|
| fail <= 2 gates | 4/4 | 4 |
| occupancy must pass | 4/4 | 4 |
| **breadth AND occupancy must pass** | **4/4** | **4** |
| breadth must pass | 4/4 | 6 |
| depth_cv must pass | 3/4 | 3 |
| identity must pass | 2/4 | 6 |

Breadth and occupancy chosen over a bare count because they mean something: a finding that
does not cover its reference is not a detection. `--actionable-requires`, default
`breadth,occupancy`, empty restores the old behaviour. NOT applied to control rows - a
control is never actionable and its content is QC, reclassified to `control_background`;
applying the rule there lost that QC signal while preventing no false positive, which
`tests/test_pipeline.py tier:NTC` caught.

Net, with the credibility gate off:

| batch | before | after |
|---|---|---|
| zymo_rl100 | 3 TP / 12 FP | **3 TP / 2 FP** |
| zymo_rl150 | 1 TP / 0 FP | 1 TP / 0 FP |
| lod_v4 (simulated) | 6 TP / 0 FP | **unchanged** |
| d12fine (simulated) | 9 TP / 0 FP | **unchanged** |

### What they exposed and nobody has fixed

  1. ~~**Two false positives survive, both *L. monocytogenes***~~ **— NOT false positives.
     Investigated 2026-10-08 and the truth model was wrong, not the pipeline.** Those two
     calls are *L. monocytogenes* at breadth 0.74 in the *B. subtilis* library and 0.21 in
     *S. cerevisiae*, and Kraken2 says both libraries genuinely contain it: **0.28%**
     (59,233 reads) and **0.06%** (5,994 reads). The alignments sit at 96.5% identity over
     101 bp, far too high for cross-genus mis-mapping - *B. subtilis* to *Listeria* would be
     ~75-80% in conserved regions at best - and reads land on four *Listeria* species at
     once. It is real DNA: index hopping or lab carry-over, which single-organism libraries
     sequenced alongside community standards are expected to carry.

     The breadths match coverage theory for those read counts, which is the check that
     settles it: 59,233 x 101 bp over a 2.9 Mb genome is ~2.2x, predicting breadth ~0.89
     against 0.74 observed; 5,994 reads is ~0.21x, predicting ~0.19 against 0.21 observed.

     **So breadth is discriminating real trace contamination from cross-mapping noise, which
     is exactly what D18's tier-ladder rule asks it to do.** The calls the rule REJECTED have
     no Kraken2 support at all - *B. anthracis* in the *S. aureus* library is 134 reads
     (0.00%) at breadth 0.0012, *T. gondii* in *S. cerevisiae* is 408 reads at 0.0001. Real
     contamination spreads over a genome; noise piles on conserved loci.

     The corrected accounting for `zymo_rl100` after the fix is therefore **3 true positives,
     2 correct detections of trace contamination, and 0 false positives** - where before the
     fix the same data gave 10 additional calls with no Kraken2 support behind any of them.
     The earlier figures in this entry and in docs/validation.md counted those two
     contamination detections as errors; they are not.
  2. **Near-identical pairs are silently attributed, not flagged.** *E. coli* (3,349,860
     reads, breadth 0.667) went to *Shigella* at ratio 0.9032, and *C. neoformans*
     (5,870,155 reads, breadth 0.642) to taxid 40410 at 0.7072 - organisms present at 100%
     abundance, reported as absent. The pack's own notes say *Shigella* is genomically
     *E. coli* and species-level separation from metagenomes is unreliable. The honest
     outcome is "cannot separate", not a silent attribution to one side.

     **Half-fixed: the neighbour is now NAMED, and that is all.** The finding read
     "reads favour neighbour taxid 620" - a number, nowhere resolved to a name in any
     artifact an analyst reads. `score_evidence.py` now loads the pack's refs manifest
     (`--refs-manifest`, passed by the pipeline; absent degrades to the bare taxid and
     changes no tier) and the same finding reads "reads favour Shigella [taxid 620]
     (near-neighbour reference) (ratio 0.9032 ...)", with `best_neighbour_name` as its
     own column. Re-scored over `runs/zymo_rl100`, 45/45 attributed rows carry a name.

     THE DECISION ITSELF IS UNCHANGED AND STILL WRONG FOR THESE TWO PAIRS. Naming
     makes the call auditable; it does not make it right. The measurement says these
     are not D12 pedestals - 620 covers 0.7328 of its reference and 40410 covers
     0.8607, so competitive alignment is working exactly as designed, and the reads
     genuinely do fit the neighbour better. They are taxonomic siblings: *Shigella* is
     a pathotype within *E. coli*, and *C. deneoformans* was *C. neoformans* var.
     *neoformans* until 2015. Whether a sibling-complex attribution should count as a
     detection of the watchlist taxon is a POLICY question over all 89 taxa, not a
     threshold to tune, and it is not answered here.
  3. **MARKER_SEARCH_AA had no time limit of its own** and inherited the 4 h base default
     while being the most expensive process in the pipeline: 1 m 41 s - 8 m 13 s on
     simulated Illumina, 33 m - 2 h 42 m on real Illumina, and killed at 4 h on ONT, whose
     4 kb reads make a six-frame protein search far larger. Now `24.h * task.attempt`.
  4. **The retry that should have covered that could not fire.** `errorStrategy` lists 143
     so an over-time task retries with more time, but Nextflow reports "process hasn't
     exited" and the trace records `status=RUNNING, exit='-'` - so `task.exitStatus` is
     NULL, `null in [...]` is false, and the run finished instead of retrying. A null status
     now retries, bounded by `maxRetries`.

**`zymo_ont` has not been run to completion.** It failed on (3) and (4) above and was
deliberately not restarted; the fixes are in but the ONT measurement does not exist yet.

**The limit that remains.** This is one standard, one lot, one study, ten organisms, four of
them on the watchlist. It is enormously better than one simulated background, and it is not
a clinical specimen: no host, no degraded DNA, no environmental diversity. Every threshold
in this pipeline is still calibrated against material that flatters it.

---

## D17 — The marker tier is permissive by omission, found 2026-10-06

Three findings from the third adversarial review. All change behaviour across many taxa, so
they are recorded rather than patched.

### D17a — 31 of 89 watchlist taxa bypass marker gating entirely, and curating a rule makes the pipeline STRICTER

`score_evidence.py` reaches `confirmed_in_silico` from the alignment gates alone, then runs
the veto / vaccine / `confirm_requires` block only `if spec and mr is not None`. When a taxon
has NO block in `marker_rules.yaml`, the whole block is skipped silently - nothing is appended
to `reasons`.

Measured on the shipped pack: 58 rule blocks against 89 watchlist entries, so **31 taxa have
no block**, 29 of them `notifiable` and 2 `emerging` - Monkeypox virus, Dengue, Zika, West
Nile, Japanese encephalitis, three orthohantaviruses, *Salmonella enterica*, *Plasmodium
vivax*, *Trypanosoma cruzi*, *Histoplasma capsulatum*, *Cryptococcus neoformans* among them.

The incentive is backwards, and the file says so about the other case: a taxon that HAS a
block but no `confirm_requires` is capped at `supported` with "no confirm_requires rule for
this taxon - cannot reach confirmed_in_silico by design". A taxon with no block at all is
confirmed on alignment alone, with no such line and no caveat. Curating a rule can only make
the pipeline stricter; omitting one makes it permissive.

Options: cap every taxon without a block at `supported` (honest, and would reclassify
confirmations across 31 notifiable taxa at once); or require a block per watchlist entry and
fail the run without one (loud, and blocks the pipeline on a curation backlog); or add an
explicit `no_markers_available: true` so the omission has to be deliberate and is visible in
`reasons`.

### D17b — the HMM tier ignores the catalogue's own thresholds

`marker_evidence.py` gates the HMM tier on `hmm_bits >= --hmm-min-bits` (default 25.0) alone.
`hmm_coverage` is computed and recorded but never gated, and the per-marker `min_coverage` /
`min_identity` the catalogue curates are not applied to this tier at all. Since a marker is
PRESENT if ANY applicable tier passes, one 25-bit domain hit makes the marker PRESENT
regardless of its curated thresholds - and 25 bits on a six-frame translation of short reads
is a weak bar. Direction of harm is false PRESENT, which fires a wrong veto or a wrong
`vaccine_if`.

Same file: the protein tier's coverage gate is skipped entirely when `aa_length == 0`, which
happens for any marker whose header in `markers.aa.faa` does not follow `marker|<id>|...`,
leaving identity as its only threshold.

### D17c — `resolution_note` is write-only

`resolve_watchlist.py` records synonym matches and RANK MISMATCHES as a note and writes it to
`watchlist.resolved.tsv`, where nothing reads it. A watchlist row whose name resolves to a
genus instead of the declared species silently expands `descendants()` over the whole genus,
so every species in it feeds that taxon's `assigned_reads`. The shipped pack is currently
benign - 13 synonym matches, no rank mismatch - but the guard rail is write-only. Note also
that `resolve_one`'s tie-break puts the rank test ahead of the scientific-name test, so a
synonym-matched node can beat the scientific-name node when their ranks differ.

---

## D16 — Two defects that need a decision rather than a patch, found 2026-10-06

Both surfaced in an adversarial review of the core scoring code. Neither is fixed, because
each changes behaviour across all 89 taxa and the right setting is not obvious from one
simulated batch.

### D16a — a deliberately contaminated NTC poisons every OTHER scenario in the same batch

`contamination_index_hop` (benchmarks/scenarios/human_core.yaml) ships with `include_ntc:
true` AND `index_hop_fraction: 0.001`, which is exactly what the scenario is for: it exists to
test the `suppressed_control` path, and its own expectation asserts the control shows the
target at `control_background`. That is correct in isolation.

It is not correct in company. `score_evidence.py` builds `control_rpm` across the WHOLE RUN,
scenario-agnostic. The NTC is 200,000 reads with no host, so ~200 hopped *Y. pestis* pairs
land at roughly **1,000 RPM** in a `role=control` row, and `control_rpm[632]` becomes 1,000
for every sample in the batch. A post-depletion sample sits near 16.7x its pre-depletion
level, so anything below about **300 RPM** pre-depletion is tiered `suppressed_control` -
including `lod_gram_negative` at 10 and 100 RPM, in a different scenario, whose *Y. pestis*
is genuinely present. The LoD would read 1000 RPM for a reason that has nothing to do with
detection, while the hop scenario's own expectation passes.

Not yet observed because every batch built so far used `--scenario` subsets that never put
the two together. A batch built from the whole file would show it.

The 2026-10-06 guard refuses `include_matrix_blank` + `index_hop_fraction` for the same
reason, but refusing `include_ntc` + `index_hop_fraction` would delete the scenario's purpose.
The options differ in what they cost:

  1. **Scope `control_rpm` by scenario.** Correct, and the coverage TSV has no scenario
     column - it would have to be threaded from the samplesheet through meta, or parsed out
     of sample names, which is the fragile option.
  2. **Exclude deliberately-hopped controls from `control_rpm`** for the hopped taxon only.
     Needs the simulator to mark them, so the pipeline would be trusting a benchmark-only
     annotation.
  3. **Refuse to build a batch that mixes a hop scenario with others**, leaving the scenario
     runnable alone. Cheapest, and narrows what a single batch can test.

### D16b — `mapped_fraction_of_assigned` mixes mates with pairs, so its floor is 0.25 on Illumina and 0.50 on ONT

KRAKEN2 runs `--paired` and METABULI `--seq-mode 2`, so stage-1 `assigned_reads` counts
PAIRS. `coverage_stats.py` counts primary records, which are MATES. The ratio is therefore
~2x inflated on Illumina and correct on ONT, so the pack's `min_mapped_fraction_of_assigned
= 0.50` is effectively 0.25 for paired data. The gate is live - 111 of 202 rows in `lod_v4`
fall below it - so this is a miscalibration by platform, not a dead threshold. RPM is
unaffected (mates over mates) and control-vs-sample comparisons stay like-for-like within a
platform. `parse_screen --min-reads 10` has the same shape: 10 pairs on Illumina, 10 reads on
ONT.

Fixing it means either halving the mate count for paired data or doubling the pair count, and
either way the 0.50 floor was CALIBRATED against the inflated figure - so correcting the units
without re-deriving the threshold would tighten the gate by 2x on every Illumina sample at
once. The units fix and the threshold re-derivation have to land together, against a batch
that can measure what the new floor rejects.

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
false confirmations). Cost, CORRECTED 2026-10-06: **2,756.3 -> 3,041.0 cpu-minutes** for 6
samples, Metabuli contributing 137.6. The figures first recorded here, "153.2 -> 173.7
CPU-minutes", were wall-clock: `evaluate_benchmark.py` summed `realtime` and never read
`%cpu`, so a 24-thread arm was charged a 24th of its cost in the table built to weigh cost
against accuracy. The conclusion is unchanged and the ratio is slightly larger than first
reported. Combined with the
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
