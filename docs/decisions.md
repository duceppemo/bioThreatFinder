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

## D25 — The Filoviridae entry, measured: it costs nothing and recovers the agent, but only after both geometry statistics were replaced

2026-10-14. D24 built the machinery; this is the entry and its price. Two batches were
run against it: `divergent_v1` (6 samples, does the lead fire) and `lod_v4` (24 samples,
9 of them `negative_host_only`, does it misfire).

**The benefit.** `Filoviridae` fires on 3 of 3 divergent replicates at 4,017 RPM and
5.43x depth, with its member covered end to end. The baseline it replaces is **zero
findings on all six samples** — complete silence on a correctly classified agent. The
descendant-aware join recovers the species: `Bundibugyo ebolavirus → LEAD via
Filoviridae`, 3 of 3, and `lead.tsv` reads `samples_led 3, covering_truth 3,
covering_nothing 0, recovered 1`. Tier PPV for `group_lead` is 1.000.

**The cost is zero.** On `lod_v4` the entry produces no findings row in any of the 24
samples and is not even nominated at stage 1 in any of them — the family node gets no
reads to roll up where no filovirus exists. `lod.tsv` and `vaccine.tsv` are
byte-identical to the pre-entry baseline. Not a cache artefact: the coverage files carry
the new `best_reference_*` columns, so `COVERAGE_STATS` re-ran under current code.

**Two defects had to be fixed first, and either would have priced the entry as a
failure.**

*The join matched nothing.* The benchmark truth says `Bundibugyo ebolavirus`. NCBI's
scientific name for 565995 is `Bundibugyo virus`; the truth string is an **equivalent
name**. The descendants writer emitted scientific names only, so the group lead would
have read as an unlisted false positive and the species it covered would have stayed a
flat miss. `NAME_CLASS_PREFERENCE` already encodes this on the resolution side, with the
comment that the official disease lists lag NCBI by years; the same list now applies to
the output side. Cost: 134,516 → 146,744 rows, 11 MB → 12 MB.

*Both geometry statistics were diluted, not one.* The first run produced the right
finding at the wrong tier: `screened_only`, the bucket for classifier noise. `breadth`
read **0.0858, which is 18,940 / 220,597** — the covered member's share of a 12-genome
bundle, and nothing about the evidence. D24 predicted geometry might let a group row
through by accident; what happened is that it kept one out by accident. `coverage_stats`
now also reports `best_reference_breadth`, `best_reference_occupancy` and names the
reference, accumulated free in a loop that already has the contig name.

**Substituting breadth alone was not enough, and this is the part worth remembering.**
The gate admitting a row to `supported` requires breadth **and** occupancy, and occupancy
is windows over the whole bundle — diluted by the identical denominator. The first fix
read correctly, survived review by eye, and changed nothing: the row failed on occupancy
and stayed in `screened_only`. The test built from the measured row is what caught it.

The substitution is for group-rank entries **only**. A species whose bundle carries three
assemblies must keep being judged on all of them, or a target covered in one of three
reads as fully covered. Asserted in both directions.

**The published row keeps the measured bundle numbers**, so a reader comparing breadth
0.0858 against a 0.20 floor would conclude the gate is broken. A `GROUP GEOMETRY:` reason
now names the reference used and what it scored.

**Seven observations moved in `lod_v4`, and none of them are this work.** Six *Y. pestis*
and one *P. vivax*, `screened_only → attributed_to_neighbour`. No group-rank taxon is
among them, and the substitution executes only for group-rank entries, which had zero
rows in that batch — so it ran zero times. The same row carries identical numbers in both
runs (2,328 reads, ratio 0.8553, uniformity 0.1524) and differs only in verdict: the
stored `findings.tsv` says *"NOT CREDIBLE: uniformity 0.1522 < 0.5, attribution
withheld"* while this run says *"reads favour Yersinia enterocolitica (ratio 0.8553 <
1.5)"*. That is the `min_neighbour_uniformity` discrepancy recorded in D21 — argparse
defaulted to 0.50 while `nextflow.config` ships 0 — and the stored findings predate the
correction. The new run follows the shipped config.

**Regenerating the resolved watchlist also refreshed `n_descendants` for the other 89
entries**, which had been resolved against an older taxdump: *B. anthracis* 51 → 56,
*Y. pestis* 287 → 291, *F. tularensis* 117 → 123. The pipeline always resolves at runtime,
so this affects only the committed copy that the tests and `build_pack_refs.py` read.

**Operational cost paid once.** `bin/` is on the Nextflow task PATH — not staged, not
hashed. Editing `coverage_stats.py` while `lod_v4` was in flight meant a live task picked
up a half-applied edit and died with `too many values to unpack`, taking the batch with
it. The same property means `-resume` does not invalidate on a `bin/` change, which is
the quieter half of the hazard.

**Still open.** The benchmark expectation
`divergent_bundibugyo_filovirus / Zaire ebolavirus / min_tier >= screened_only` still
FAILs with *no finding*. It was written assuming the divergent reads would cross-map onto
Zaire; they do not, and all 3,992 go to Bundibugyo. The expectation encodes the old theory
of how this should be caught and should probably become an expectation on the lead — but
that is editing the benchmark's own ground truth and has not been done. The other half of
D15, the BoNT detection failure with all three markers ABSENT and every tier empty, is
untouched by any of this.

---

## D24 — A group-rank watchlist entry needs its own tier and its own join, or it is wrong in both directions at once

2026-10-14. D15's divergent filovirus is classified **correctly** and reported **not at
all**: Kraken2 and Metabuli each put all 3,992 spiked reads on *Bundibugyo virus* and roll
them up through Filoviridae, while every species-keyed watchlist row stays silent. The
obvious repair is a family-rank watchlist entry. Verifying that repair before building it
turned up two things that would each have made it useless.

**The rollup itself needs no code.** Checked against the published run and the real
taxdump, four of five links already hold: Filoviridae carries clade_reads 3,992 in the
Kraken2 report with 0 direct; `clade_or_sum()` takes the node's own clade_reads, so a
family row is credited all of them; `resolve_watchlist.py` resolves `Filoviridae` to taxid
11266 at rank `family` with 103 descendants, no warning and no rank-mismatch note; all
five filovirus taxids are inside that descendant set. `build_pack_refs.py` queries the
taxid's subtree with a reference-first pass, so a family entry fetches one curated
assembly per species.

**But a group entry scored on alignment geometry would be a species call.** Breadth for a
taxid divides the covered span by the summed `reference_bp` of every genome sharing it, so
a family bundle of eight filovirus genomes covered only in its Bundibugyo member lands
near 0.13 occupancy — just under the 0.15 floor, which is the tier we want, reached by
arithmetic that a ninth genome would overturn. Hence an explicit cap: an entry whose
`resolved_rank` is supra-specific is capped at a lead, by construction rather than by
accident.

**The cap must be an allowlist of group ranks, not `!= "species"`.** 18 of the 89 shipped
entries resolve to NCBI rank `no rank` and every one of them names a real organism —
Variola virus, Zaire ebolavirus, SARS-CoV-2. A denylist would have capped Variola.

**Two measurements then redirected the design.**

*The lead cannot reuse `screened_only`.* That tier is the **initial value** of the tier
chain — the bucket every coverage row starts in and falls back to — so it means "stage 2
established nothing". It holds 74, 59 and 23 rows on lod_v4, lod_v3_threearm and
zymo_rl100, against 6, 24 and 43 positives respectively. A lead routed there is
indistinguishable from the pile. The alternative considered and rejected was adding
`screened_only` to `POSITIVE_TIERS`, which would have turned 74 unresolved rows into
positive calls on a run that makes 6. So `group_lead` is its own tier, placed above
`attributed_to_neighbour` (which asserts the target is NOT here) and below
`indistinguishable_complex` (which narrows the claim to a named pair — a group of 103 is
the weaker statement).

*The evaluator could not have scored it.* `evaluate()` joins findings to benchmark truth
on the **exact taxon name**. The truth says `Bundibugyo ebolavirus`; the finding says
`Filoviridae`; the strings never meet. Left alone the group finding lands in `unlisted` as
if it were a false positive on a taxon nobody listed, **and** the species it covers still
scores a flat miss — wrong in both directions simultaneously, which is why the family
entry would have measured as a regression. Fixed with a descendant-aware join:
`resolve_watchlist.py --descendants-out` now writes `descendant_name` beside each taxid
(the manifest's truth columns carry names only), and `load_covers()` feeds the join. The
extra column is backward compatible — `parse_screen.py` reads that file positionally with
`len(parts) >= 3`.

**A lead is neither a hit nor a miss, and both would be a lie.** Scoring it POSITIVE would
credit a species call the pipeline explicitly declined to make, and would let one entry
collect a true positive for any of 103 descendants. Scoring it NEGATIVE erases the fact
that the right group was put in front of the analyst. So `LEAD` is its own outcome, and
**sensitivity does not improve**: a lead-covered species stays in the denominator and
stays a non-positive, with `n_lead_only` beside it. Specificity is untouched in the other
direction — `LEAD` is deliberately absent from `PRESENCE_ASSERTING`, because a lead names
no member and so accuses none of them. `lead.tsv` reports both sides: species recovered
only as a lead (each one arguably a species entry the watchlist is missing, so a number
that stays high is a gap, not a win) and leads that fired where the truth lists no member
of the group (the cost of carrying the entry).

**Ordering bug caught by its own test.** The cap first ran before the control
reclassification, which keys on `ASSERTS_PRESENCE`. A capped control row was no longer in
that set, so a contaminated control stood as an open lead in its own report instead of
being named as background. The cap now runs after it.

**Nothing in the shipped pack is above species rank, so no finding changes today**, and a
test asserts that — it is the notice if it stops being true. Re-running the evaluator over
lod_v4 (24 samples, 202 findings) against HEAD on identical inputs: `specificity`,
`tier_ppv`, `lod` and `vaccine` byte-identical, `sensitivity` and `observations` differing
only by the new columns.

**Still open, and the reason D15 is not closed here.** This is the machinery, not the
entry. Whether to actually add Filoviridae is a specificity question that needs a run of
`negative_host_only` and `divergent_v1` to price, and the half of D15 that is a detection
failure — all three BoNT markers ABSENT with every tier empty — is untouched by any of
this.

---

## D21 — D16a closed by refusing the batch, not by scoping the scorer; and control subtraction stops being invisible

2026-10-13. Fixing D16a meant first discovering that **the behaviour D16a calls a bug is
the behaviour D14 depends on**, and that the defect nobody had looked at was larger than
the one on the list.

### The fix I started building would have reinstated D14

The obvious remedy is D16a's own option 1: scope `control_rpm` by scenario so a
deliberately contaminated NTC stops suppressing other scenarios. I designed a general
version of it — an optional `control_group` samplesheet column — before finding this
comment in `negative_host_only`'s `controls:` block:

> ONE PER BATCH IS ENOUGH, and it belongs here rather than on every scenario:
> `score_evidence.py` builds `control_rpm` across the whole RUN, not per scenario, so
> this blank is what every sample in the batch is compared against.

`[MEASURED]` The comment is behaviour, not intent. Re-scoring `lod_v4` with and without
its controls — the controls are all `negative_host_only_*`, a **different scenario** from
the 15 `lod_gram_negative` samples they act on:

| variant | confirmed | suppressed | screened |
|---|---|---|---|
| all controls (today) | 6 | 36 | 71 |
| NTCs only | 6 | 18 | 81 |
| **no controls** | **24** | **0** | 69 |

Removing the cross-scenario reach moves 30 rows and hands back 15 false
`confirmed_in_silico` calls on *E. coli* plus 15 *T. gondii* false positives. **The
run-global default is load-bearing and stays.** One nuance that revises D14's framing:
under the current pack the *T. gondii* artefact lands at `screened_only` without a
control rather than `supported` — the read floor and D20's alignment work already demote
it — so the live cross-scenario exposure today is *E. coli*, not *T. gondii*.

### So D16a is closed by refusing the batch

Exactly **one** scenario in the repo declares `index_hop_fraction`, and
`build_matrix_blank_plan` already refuses `include_matrix_blank` + `index_hop_fraction`
for this exact reason. The guard is extended to the batch level: a scenario with
`index_hop_fraction` cannot share a batch with another scenario, and the error names the
`--scenario` escape hatch. It costs the ability to run the hop scenario alongside others
in one batch, which per D16a's own text **nobody has ever done**, and it cannot regress
D14 because it does not touch scoring.

**"Closed" means the benchmark manifestation, not the defect.** A real laboratory with a
genuinely contaminated NTC reproduces D16a with no benchmark anywhere near it, and
nothing here prevents that — run-global subtraction is *correct* for that lab, and what
was missing was any way to see it happening. The audit below is that remedy. Refusing
the batch closes the way this project can build the defect; it does not close the way the
world can present it. Those are different claims and this file should not blur them.

Three routes also get past a build-time guard even inside the benchmark: hand-merging two
samplesheets, running the pipeline over two separately-built batches, and re-scoring two
batches' published coverage TSVs together — the last of which is how most measurement in
this project is done. `evaluate_benchmark.py` therefore warns when a manifest contains
any `index_hop=` row **and** more than one `scenario_id`. A warning rather than a
refusal: by the time a manifest is read the reads exist and the pipeline has run, so
exiting destroys a completed evaluation without preventing anything.

`control_group` is **not** shipped. It is the better long-term abstraction, but its
trigger rule was wrong twice in two days (scenario id → `include_ntc` →
`index_hop_fraction`), the samplesheet has no schema documentation anywhere in the repo
to document it in, and it would arrive with exactly one writer and no non-benchmark
consumer. A general mechanism with one benchmark call site is a benchmark annotation
wearing a disguise. Revisit it with the `control_type` work, and reconcile the two then.

### The larger defect: control subtraction has been invisible

`[MEASURED]` Across six published batches control subtraction moved **108 rows**, and
every number in this entry had to be recovered by regexing the prose of the `reasons`
column — the artifact that breaks silently the moment a message is reworded. Three things
were not visible anywhere:

  1. **Which control decided a suppression**, in any machine-readable form.
  2. **11 suppressed rows carry MORE of the organism than the control does** — 8 in
     `lod_v4`, 3 in `d12fine`, at sample/control ratios of 1.01–1.03. The contamination
     model did not decide these. `control_ratio 5.0` did, and nothing in this file
     derives that 5.0 from anything.
  3. **Controls carry watchlist pathogens that suppress nothing yet.** `lod_v4`'s
     controls hold *Bacillus anthracis* at 40.09 RPM — a standing 200.46 RPM floor on a
     select agent — plus *N. meningitidis*, *Shigella dysenteriae*, *L. monocytogenes*
     and *S. enterica*. None suppresses anything today, because every row for them was
     taken by the read floor or by `attributed_to_neighbour` before suppression could be
     reached. A control holding a select agent and firing zero times is not a clean
     control; it is an armed one, and the pipeline said nothing.

Landed: a stderr audit table (one line per taxon ANY control raised, with the implied
floor, and `LATENT` where nothing was suppressed), a `control_subtraction` column in
`findings.tsv` and the summary TSV (`cleared:` / `suppressed:` / `exceeds_control:`, each
naming the deciding control), an explicit reason string on the exceed-control rows, and a
"Control subtraction" section in the HTML report. `[MEASURED]` Behaviour-neutral:
**0 tier changes across all six batches.**

### An argparse default that disagreed with what ships

`[MEASURED]` D18 set `min_neighbour_uniformity = 0` in `nextflow.config`. The argparse
default in `score_evidence.py` stayed **0.50**, so every manual re-score that forgot the
flag silently reproduced **pre-D18** behaviour. This is not hypothetical: it happened
during this review. Injecting one 1,000 RPM control row into `zymo_rl100` produces **5**
tier changes under the shipped configuration and **7** under the argparse default, and
the two extra rows are *E. faecalis* and *S. aureus*, which the shipped settings hand to
`attributed_to_neighbour` before suppression can see them. The 7 was reported as a
correction to this file's D16a note; **the note was right and no correction was needed.**
The default is now 0, and `tests/test_pipeline.py defaults:` asserts that every threshold
present on both sides agrees — including a check that the list of compared pairs is still
the complete list, since an uncompared new pair is how this one survived.

`[MEASURED]` A consequence worth stating plainly: **the published findings under `runs/`
are no longer reproducible at the shipped configuration.** The config moved to 0 at
`c49beb1` (2026-10-08) and `zymo_rl100` was produced on 2026-10-07 — its published
findings carry 16 `attributed_to_neighbour` rows, which a re-score reproduces at 0.50 and
turns into 45 at 0. D20's alignment floor moved too. **A re-score of a published batch is
a different measurement from the batch**, and every comparison in this file that
re-scores must state which configuration it used. The `defaults:` test removes one cause
of that drift, not all of them.

### What this did not fix

`control_rpm` is still `max()` over controls with **no prevalence requirement**, so a
single contaminated control decides — the literature wants prevalence across ≥5 controls,
flag-not-delete, which is roadmap item C. Grouping would have made this *worse*, not
better: each group holds fewer controls, usually one, so `max()` becomes "this one
control decides" by construction. `control_ratio 5.0` remains underived; the audit
table's `floor RPM` column is what would let it be re-derived. And the mixed batch that
would demonstrate D16a as a live defect rather than an injected one still does not
exist — the guard now refuses to build it, which is the point, but it means the fix is
pinned by fixtures and a synthetic injection rather than by a run.

**What would reverse it:** a real deployment with several extraction batches under one
samplesheet, where run-global subtraction is wrong for a reason that is not a benchmark
artefact. That is when `control_group` should land — with a samplesheet schema document,
a non-simulator fixture, and reconciled against `control_type`.

---

## D23 — The marker tier's false positives were a protein-identity discount, not the HMM tier; and half of D15 is not a reporting bug at all

2026-10-14. Work on D15 (`marker_lead`) and D17b (the HMM tier) measured three things
that redirect both.

### D17b's own worked example survives D17b's own fix

D17b says the HMM tier gates on `hmm_bits >= 25` alone and never applies the catalogue's
`min_coverage`. True. `[MEASURED]` Applying it changes nothing — on every HMM-PRESENT
marker in the filovirus sample, observed coverage already clears the catalogue floor:

| marker | catalogue min_coverage | observed hmm_coverage | killed? |
|---|---|---|---|
| marburg_L | 0.50 | 0.7160 | **no** |
| ebola_GP | 0.60 | 0.9972 | no |
| ebolaS_GP | 0.60 | 0.7862 | no |
| ecoli_stx2 | 0.80 | 1.0000 | **no** |

The fix the entry points at is a no-op on the case the entry cites. And no bits
threshold separates them either: `marburg_L` at 93.7 bits is *correct cross-genus
filovirus homology* on a Bundibugyo sample — verified by the negative control, where
every filovirus marker is ABSENT with empty HMM fields. `marburg_L` carries the lowest
floors in the catalogue precisely because its note reads "most conserved filovirus gene;
best HMM target". **Any threshold tight enough to suppress it on Bundibugyo also
suppresses it on a divergent Marburg lineage detected only on L, which is the case it
was curated for.** The HMM tier is not broken; the reporting unit is — see below.

### The real false-positive route is the protein tier, and the catalogue already had the answer

`aa_min_id = max(0.50, min_id_catalogue - 0.15)` subtracted 15 identity points from every
marker regardless of curation. `[MEASURED]` On
`divergent_bundibugyo_filovirus_illumina_r1` — one virus, no other watchlist pathogen —
that made **six** markers PRESENT off unrelated background, every one a conserved
housekeeping gene:

| marker | catalogue min_id | old effective floor | observed aa_identity |
|---|---|---|---|
| mtb_rpoB | 0.98 | 0.83 | 0.877 |
| cdiph_rpoB | 0.97 | 0.82 | 0.900 |
| bholm_recA | 0.97 | 0.82 | 0.860 |
| mtb_katG_inhA | 0.97 | 0.82 | 0.840 |
| cburn_icd | 0.95 | 0.80 | 0.877 |
| rprow_gltA | 0.95 | 0.80 | 0.833 |

The premise ("protein floors are lower than nucleotide by nature") is true and the
implementation inverted the curation that already expressed it: cross-matching
housekeeping genes are floored at 0.95–0.98, while the divergent-virus block is floored
at 0.78–0.90 **on purpose**, with a note saying identity floors there are deliberately
low because the HMM tier carries detection. Honouring the catalogue removes all six and
keeps both real filovirus hits — `ebola_NP` 0.880 and `ebolaS_NP` 1.000 against a
deliberate 0.85 floor.

`[MEASURED]` Corpus-wide, across **143 marker_states files**: **828 PRESENT observations
lose their only passing tier, spanning 13 markers.** Six housekeeping genes account for
820 of them. The other 8 are all false positives on real Zymo data, each with
`nt_reads = 0` — including `banth_pxo2_capA`/`capB` in a mock community containing no
*B. anthracis*, and `ypest_ppcp1_pla` in a *Salmonella* isolate. **Zero true positives
are lost.** One of the six, `mtb_katG_inhA`, is `report_always`, i.e. it was eligible to
surface as "isoniazid resistance" on samples with no mycobacteria.

The slack is now a declared argument defaulting to **0**, with a separate, smaller ONT
allowance that does not compound with the nucleotide relaxation.

### Why this is behaviour-neutral today and would not have been tomorrow

`[MEASURED]` Re-scoring `lod_v4`, `d12fine` and `zymo_rl100` with the new floor applied:
**0 tier changes and 0 marker_verdict changes.** The reason is the whole point —
*M. tuberculosis*, *C. diphtheriae*, *C. burnetii* and *R. prowazekii* have **no coverage
rows** in those batches, and marker evidence can only gate a finding that already exists.
The 828 false PRESENTs are inert **only because markers cannot raise findings**. D15
proposes exactly that change. So the tier fix must precede any `marker_lead`, and the
ordering in the backlog ("D15 must land with D17b") was right about the dependency and
wrong about which defect.

### Half of D15 is not a reporting failure

`[MEASURED]` On `divergent_clostridium_bont_illumina_r1`, **all three BoNT markers are
ABSENT** — `cbot_bont`, `cbot_ntnh`, `cbut_bont`, each searched via `dna,hmm,protein`,
every tier empty. There is no marker evidence to raise under any rule. D15's second
scenario is a **detection** failure, not a reporting one, and needs its own
investigation starting with `cbut_bont`'s `partial` manifest status.

### `marker_class` is not a usable gate, which kills the obvious `marker_lead` rule

The tempting rule was "raise a lead only on a PRESENT `discriminating` marker" — clean on
both samples tested. `[MEASURED]` It does not generalise:

  - **24 of 89 watchlist taxa have a `discriminating` marker at all**, so 65 could never
    lead. That is D17a's silent-coverage-gap shape with the sign flipped.
  - **4 of those 24 have only a MULTI-COPY one** — *B. pertussis* (`IS481`),
    *C. burnetii* (`IS1111`), *B. melitensis* (`IS711`), *B. mandrillaris*. The
    catalogue header forbids exactly this: "multi markers must never be used alone", and
    `IS481`'s own note calls it "the classic false positive".
  - The classes are not sorted on cross-match risk. `required` holds both housekeeping
    genes (`rpoB`, `recA`, `gltA`) and conserved viral cores (`marburg_L`, `ebola_NP`);
    `discriminating` holds surface antigens and insertion sequences.

### What is NOT being changed, and the open decision

`marker_lead` is **not** built. For the filovirus sample the true answer is "a divergent
filovirus is present"; `marker_lead` would report **Zaire ebolavirus** and **Sudan
ebolavirus**, two species the scenario's own `truth.absent` names. Two wrong species is
actionable in the wrong direction in a way silence is not. The finding needs to say
"filovirus, species undetermined", and it cannot, because the resolved watchlist is 71
`species` + 18 `no rank` with nothing above species — **which is the gap D15 originally
named and the re-diagnosis discarded.** That decision (a lead-tier-only family rank entry
versus closing D15 won't-fix) is open and is recorded as open rather than pre-empted by
shipping a tier.

Also landed: `hmm_qlen` and `aa_length` are emitted into `marker_states.tsv`. Both were
computed and discarded, so neither `hmm_bits` nor `aa_coverage` could be interpreted from
the artifact — a bitscore is roughly linear in profile length, and `marburg_L`'s ~2,100 aa
profile against `ebola_GP`'s ~676 aa is the only thing distinguishing 94 bits from 106.
And `marker_rules.py --validate --watchlist` now cross-checks the three hand-maintained
files that agree only on literal strings; it reports the two real asymmetries today
(*Brucella suis* has a rule block but no marker of its own, so a PRESENT call on
*B. melitensis*'s markers raises it too; *Monkeypox virus* has a marker and no block).

**What would reverse the protein change:** a true positive demonstrably lost to a
catalogue floor. The ONT path is where to look first — translated identity is less
error-sensitive than nucleotide, but no ONT batch has ever been scored.

---

## D22 — `control_ratio = 5.0` is underived, and now there is an instrument that could derive it

2026-10-13, split out of D21 because it is a threshold with no provenance rather than a
defect with a fix.

A finding is tiered `suppressed_control` when `sample_rpm < control_ratio x control_rpm`,
with `control_ratio` defaulting to **5.0**. `[REASONED]` Nothing in this file, in
`nextflow.config`'s comments, or in the pack derives that 5. It is not traceable to a
measurement, a guideline or a citation; it is a round number that has been load-bearing
since the first run.

`[MEASURED]` What it decides, from the audit added in D21:

  - **11 suppressed rows across `lod_v4` and `d12fine` carry MORE of the organism than
    the control does** — sample/control ratios of 1.01–1.03. The contamination model did
    not decide these. The margin did, alone.
  - The margin sets an implied floor per taxon per run that nobody chose: *E. coli* at
    **1,665,171 RPM** in `d12fine`, *T. gondii* at 21,204, *B. anthracis* at **159** —
    the last on a select agent, from a control holding 31.83 RPM.
  - A 5x margin means a control at 1,000 RPM owns everything below 5,000 RPM. All 9 of
    `d12fine`'s genuine *Y. pestis* detections sit at 2,149–2,925 RPM, inside that band.

**Why it is not being changed here.** Lowering it admits contamination; raising it
suppresses more true positives; and this entry cannot say which, because the project has
never measured the distribution of sample/control ratios for *known* contamination
against *known* presence. Changing an underived number to a different underived number
is not progress.

**What would derive it.** The audit table's `floor RPM` column is the missing instrument:
it makes every implied floor visible per run, so a batch containing both genuine
low-level positives and genuine carry-over contamination would show where the two
populations separate — the same shape of argument D12 used for neighbour uniformity, and
the same caution, since that separation looked clean on five bundles and narrowed to a
0.07 gap across 182 rows. Until then 5.0 stays, declared rather than justified.

**What would reverse it:** any measurement of the ratio distribution; or a single
documented case of a true positive suppressed by the margin alone, which the
`exceeds_control` marker will now surface.

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
min_read_identity` silently dead. ✓ **D16a**, closed by refusing the batch rather than
scoping the scorer (D21) — the scoping fix was measured to reinstate D14. **D17a**, which is dormant only because nothing can
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

### D17b — the HMM tier ignores the catalogue's own thresholds ⚠️ *true, and the fix it implies is a no-op; the real route was the protein tier — see D23*

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

**2026-10-14, D23.** Measured, both halves of this entry are weaker than they read. Applying
the catalogue `min_coverage` to `hmm_coverage` kills **nothing** — `marburg_L` observes
0.7160 against a 0.50 floor, `ecoli_stx2` observes 1.0000 against 0.80 — and `marburg_L`'s
94 bits turns out to be correct cross-genus filovirus homology rather than noise, so no bits
threshold separates it from `ebola_GP` without also suppressing the divergent-Marburg case
the marker exists for. The `aa_length == 0` bypass is real but **latent**: every row with
`aa_hits > 0` also has `aa_coverage > 0` across the corpus. The live defect was the protein
tier's undeclared blanket `-0.15` identity discount, which D23 fixes and which accounts for
828 false PRESENT observations across 143 files.

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

### D16a — a deliberately contaminated NTC poisons every OTHER scenario in the same batch ✅ *closed 2026-10-13 by option 3; see D21*

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

**Option 3 was chosen (D21).** Option 1 was designed, built as far as the design, and
abandoned on measurement: scoping `control_rpm` away from the whole run reinstates D14,
because `negative_host_only`'s matrix blank deliberately relies on the run-global reach —
re-scoring `lod_v4` without its controls returns 30 rows, 15 of them as false
`confirmed_in_silico`. The three options above are therefore not equal-cost alternatives:
option 1 trades this defect for a worse one. The batch-level guard lives next to the
`include_matrix_blank` + `index_hop_fraction` refusal it mirrors.

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

## D15 — A divergent agent is invisible because it is classified CORRECTLY, and the watchlist registers its sibling as a suppressor ⚠️ *the divergent-filovirus half is CLOSED: see D24 for the lead tier and descendant-aware join, D25 for the measured entry (3/3 recovered, zero specificity cost). The BoNT detection failure is still open.* *the amino-acid arm does not fix this, and was tested on the scenario built to favour it*

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
  **2026-10-14 (D23): option 3 is measured and does not work as stated; the entry's
  original diagnosis — no reporting unit above species — is back as the live question.**
  The filovirus markers DO fire (`ebola_NP`, `ebola_GP`, `ebolaS_NP`, `ebolaS_GP` all
  PRESENT, `confirm_requires` satisfied for both Zaire and Sudan) and are discarded
  because `score_evidence.py` builds findings only from coverage rows. But a
  `marker_lead` tier would report **Zaire** and **Sudan ebolavirus**, two species this
  scenario's own `truth.absent` names, for a specimen containing neither. And it does
  nothing for the second scenario: all three BoNT markers are **ABSENT**, so that half
  is a detection failure, not a reporting one. Option 1 is the open decision.

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

---

## D26 — A tier that was never built was being reported as searched, on 13 markers, almost all of them the divergence-catching tier

2026-10-14. Found while diagnosing the BoNT half of D15. It is not what made
BoNT fail, but it is the reason that failure was hard to read.

**The defect.** `tiers_supplied` answers "did this run receive an input file for
the tier". It was being used as if it answered "does this marker have anything in
that file", and `tried` was appended unconditionally inside each tier's branch —
including the `else` arm where the lookup missed. A marker declaring
`dna+protein+hmm` with no HMM model therefore emitted

    tier_used = dna,hmm,protein
    reason    = "searched via dna,hmm,protein; no tier met thresholds"

The row asserted that the most sensitive tier looked and found nothing, when the
model did not exist.

**Scope on the shipped pack: 13 markers.** Eleven are missing the HMM tier —
`cbut_bont`, `ebola_NP`, `ebolaS_NP`, `marburg_NP`, `marburg_GP`, `junin_GPC`,
`rvfv_L`, `rvfv_Gn`, `nipah_L`, `nipah_G`, `eeev_nsP`. Every one is a
divergence-catching marker on a high-consequence agent, which is exactly the set
`markers.tsv` argues the HMM tier exists for: *"HMM tier is essential — novel/mosaic
subtypes are regularly described and would be missed by exact matching."* Two —
`orthopox_ATI` (Variola) and `mpox_C3L` (Monkeypox) — are missing the **nucleotide**
tier instead, which is the worse direction: a Variola marker whose reason claimed a
nucleotide search that never happened.

**The fix.** The manifest is authoritative per tier as well as per marker, via
`n_nt_seqs`, `n_aa_seqs` and `has_hmm`. An unbuilt tier is dropped from
`evaluable`, excluded from `tier_used`, named in a new `tiers_not_built` column
and appended to the reason as `NOT SEARCHED - declared but not built for this
marker`, on PRESENT rows as well as ABSENT ones. A per-run stderr block lists
every gapped marker once.

**One verdict moves, deliberately.** A marker whose *only* declared tier is
unbuilt was ABSENT — nothing was searched, and that silence was being reported as
evidence against the organism. It is now UNKNOWN, which can neither promote nor
veto. This is the same class of error the module was written to prevent and
already records for `mtb_RD1`: an untested marker reading as a clean negative and
turning a real TB detection into a BCG stand-down.

**Everything else is unchanged, and that is asserted rather than assumed.** Re-run
against the committed `divergent_clostridium_bont_illumina_r1` state file: 163 rows
before and after, **zero state changes**, exactly the 13 `tier_used` strings
corrected. An unbuilt tier has no hits either way, so a fix that moved a verdict
would have been changing something else as well.

**Backward compatibility is load-bearing.** A manifest predating the three columns
must not read as all-gaps and silently turn the whole catalogue UNKNOWN. Column
absence from the *header* means "not recorded, assume built"; a column present and
empty means no. Pinned in both directions.

Eleven new `tiergap` checks. `curation_required.tsv` does not list any of the 13,
which is a separate gap in `build_marker_db.py` and is not addressed here.

---

## D27 — The curation to-do list omitted 18 of the 31 markers that needed curating, including every retrieval failure

2026-10-14. Follow-on from D26. The gaps that fix made visible per-run should
never have needed finding by hand; `curation_required.tsv` exists to be that
list, and it was filtering out most of its own subject.

**The filter was `status in (curation_required, motif_todo)`**, which dropped two
categories.

**`failed` — never listed at all, 5 markers.** A retrieval failure evaluates
UNKNOWN downstream exactly like a `curate_required` row, which is the stated
reason the file exists; it was simply not in the tuple. Three of the five are
load-bearing in `marker_rules.yaml`:

| marker | role in the rules |
|---|---|
| `bthai_ara` | `veto_if` for **both** B. pseudomallei and B. mallei |
| `henipa_cedar` | `veto_if` for **both** Nipah and Hendra |
| `bruc_S19_eri` | half of the Brucella `vaccine_if` |

A veto that can never fire is the quiet kind of gap. UNKNOWN makes the rule
indeterminate, which holds a finding at `supported` rather than confirming it —
the safe direction, so nothing crashes and no call looks wrong. It costs
specificity silently: a B. thailandensis or Cedar-virus cross-match cannot be
vetoed down, because the marker that would do it was never retrieved.

**`partial` — excluded deliberately, 13 markers.** The reasoning was that these
rows do evaluate, on the tiers they have. They do, and on fewer tiers than
`markers.tsv` declares, with the missing tier recorded nowhere but a manifest
`note`. Eleven are missing the HMM tier — `cbut_bont`, `ebola_NP`, `ebolaS_NP`,
`marburg_NP`, `marburg_GP`, `junin_GPC`, `rvfv_L`, `rvfv_Gn`, `nipah_L`,
`nipah_G`, `eeev_nsP` — the tier the catalogue says exists to catch a divergent
allele. `orthopox_ATI` (Variola) and `mpox_C3L` (Monkeypox) are missing the
nucleotide tier instead.

**The fix.** `missing_tiers()` derives the gap from `evidence_type` against
`n_nt_seqs` / `n_aa_seqs` / `has_hmm`, mirroring D26's rule in
`marker_evidence.py` exactly — if the two drift, the build's to-do list and the
run's reported gaps disagree about which markers are instruments. A row is listed
if it is untestable **or** has a gap; a new `missing_tiers` column says which; a
`DEGRADED` line per marker goes to the build log, which is what an operator
actually reads.

**The header comment was also wrong and is rewritten.** It said *"Every marker
listed here evaluates as UNKNOWN downstream"*, which is false for a `partial`
row. The two kinds are now distinguished in the file itself.

13 → 31 rows on the shipped pack. Five self-test checks, including that an
undeclared tier which happens to exist is **not** a gap — `has_hmm: yes` on a
`dna+protein` marker is spare capacity, not a missing instrument.

The 18 underlying gaps are **not fixed** here; they are now tracked. `cbut_bont`
names its own obstacle — *"only 2 protein sequences"*, below the builder's
alignment floor — so the HMM set needs more sequence, not a re-run.

---

## D28 — A second expectation axis, because a findings row only ever exists for a watchlist taxon

2026-10-14. The gap D15's BoNT half was actually sitting on, found while
diagnosing it: the scenario could not state its own claim.

**What the benchmark could assert, and what it could not.** Every expectation
in `human_core.yaml` is `taxon:` plus `min_tier` / `max_tier` /
`min_tier_at_rpm` — a claim about a row in `findings.tsv`. That table only ever
carries watchlist taxa. `divergent_clostridium_bont` spikes *C. butyricum*,
which is a `near_neighbour` and deliberately **not** on the watchlist, so no
findings row can ever name it. Its only expressible expectation was

    - taxon: Clostridium botulinum
      max_tier: marker_veto

an upper bound. `tier_rank(None) == -1`, so "no finding" satisfies any cap: the
scenario **PASSed on total silence**, on every run, while the toxin it exists to
test went unreported. Its own description says *"Correct behaviour is to report
the toxin and name the right organism"*, and nothing checked that.

The same hole covers every toxin-in-a-neighbour scenario — Stx in a non-O157
*E. coli* is the other obvious one.

**The fix is a second axis, not a wider first one.**

    - marker: cbut_bont
      state: PRESENT        # PRESENT | ABSENT | UNKNOWN

fed by `--marker-states`, the per-sample files `marker_evidence.py` already
writes. Samples come from the **manifest**, not from `obs`: a scenario whose
spike is untargeted has no observations at all, so grading off `obs` would have
reproduced the same blindness one level down.

**Graded over every sample, all of which must agree** — the rule `min_tier`
already uses. A marker that fires on two replicates of three is not a working
marker.

**UNKNOWN is a database gap, not a result**, and is excluded from grading
exactly as `no_reference` is for tiers: an unbuilt marker SKIPs rather than
FAILs. That symmetry is the contract `marker_evidence.py` is built on — an
untested marker may not read as a clean negative, and it may not read as a
detection failure either. `state: UNKNOWN` grades normally, because then the gap
*is* the assertion, which is what makes D26's 13 markers pinnable. D26's
`tiers_not_built` is surfaced in the detail, so a SKIP explains itself.

Three failure modes are ERROR rather than FAIL, because they are authoring
mistakes that would otherwise be expectations which cannot fail — the
`check_tier_name` lesson applied to the new axis: an unknown state name, a
missing `marker`/`state` pair, and a marker id that appears in no state file at
all. A sample missing a row its siblings have is a genuine FAIL.

**Measured against the real run.** With the expectation added to a scratch copy
of the scenario file and the committed `divergent_v1` marker states:

    divergent_clostridium_bont  Clostridium botulinum  max_tier                 <= marker_veto  no finding  PASS
    divergent_clostridium_bont  Clostridium butyricum  marker_state[cbut_bont]  PRESENT         ABSENT      FAIL
    divergent_clostridium_bont                         marker_state[cbot_bont]  ABSENT          ABSENT      PASS

The scenario goes red, correctly, naming all three replicates — and the old
taxon expectation is untouched, which is asserted rather than assumed.

**The scenario file is NOT edited here.** The machinery exists; adding the
expectation is a change to benchmark ground truth and is the user's call, as is
the fixture strain. Both are recorded in the BoNT diagnosis: the spike is
*C. butyricum* DSM 10702, the non-toxigenic type strain, which carries no BoNT
gene (zero minimap2 hits, zero tblastn hits at 1e-3, against a validated
positive control). Adding the expectation without fixing the fixture makes the
scenario correctly red; fixing the fixture without the expectation leaves it
unable to fail.

Ten self-test checks, mutation-verified, plus ten CLI-level `markerexp` checks
covering the flag, the file reader and the exit code. `run_bakeoff.sh` passes
marker states from the same arm the reference findings came from.

---

## D29 — The BoNT scenario now states its claim, and is correctly red

2026-10-14. D28 built the axis; this is the expectation, added on the user's
explicit direction since it changes benchmark ground truth.

```yaml
- marker: cbut_bont
  state: PRESENT
  taxon: Clostridium butyricum
```

**Only `cbut_bont` is asserted, and the omissions are the substance of this
entry.** A first draft also expected `cbot_bont: ABSENT`, which would have been
wrong as ground truth. BoNT/E is shared between the two species, `cbot_bont`
lists `botE` among its aliases, and `marker_rules.yaml` is built around exactly
that overlap: `veto_if: cbut_bont` exists so a `cbot_bont` hit can be
*attributed* to the right organism rather than suppressed. Asserting it absent
would contradict the rule the scenario exists to exercise. `cbot_ntnh` is the
same case — NTNH is co-located in the toxin cluster in both species.

**`min_tier: marker_veto` on *C. botulinum* is NOT added.** The description says
the veto should name the right organism, but a veto only produces a row if the
classifier nominates *C. botulinum* at all, and there is no evidence from any
run that it does. That expectation would be a guess about classifier behaviour
dressed as ground truth.

**Measured against the committed file and the real `divergent_v1` run:**

    divergent_clostridium_bont  Clostridium botulinum  max_tier                 <= marker_veto  no finding  PASS
    divergent_clostridium_bont  Clostridium butyricum  marker_state[cbut_bont]  PRESENT         ABSENT      FAIL

2 pass, 2 fail, 26 skip. **The red is correct and is not a regression.** The
spike is *C. butyricum* DSM 10702, the non-toxigenic type strain, which carries
no BoNT gene at all — zero minimap2 hits for 14 marker nucleotide sequences,
zero tblastn hits for 12 marker proteins at e < 1e-3, against a positive control
(`ecoli_stx1` vs the background *E. coli*) that returned hits at e < 1e-110. The
benchmark is now reporting a real fixture gap instead of a green light. Swapping
the spike to a toxigenic strain is a separate change; the simulator already
accepts `genome: <accession|path>` on a spike.

**Two static guards added**, the same orphan check the scenario taxa already
get: every `- marker:` must name an id in `markers.tsv`, and every `state:` must
be a real state. The evaluator reports both as ERROR at run time, but the lesson
`check_tier_name` records is that a misspelled assertion silently removes the
gate it was written to be — catching it in the suite means a typo fails in
seconds rather than hours into a benchmark run.

---

## D30 — The family relaxation that was written for these exact markers had been unreachable for its whole life, and making it reachable was the wrong fix

2026-10-15. Follow-on from D27, which tracked 18 marker gaps and fixed none of
them. Eleven were the HMM tier missing on divergence-catching markers for
high-consequence agents. This entry closes nine of them, and the interesting
part is the two it deliberately leaves open and the regression it nearly
shipped.

**The dead remedy.** `load_rankedlineage()`'s docstring names its own purpose:

> A profile needs at least three sequences to mean anything, and the first full
> build produced only 1-2 for ebola_NP, ebolaS_NP, marburg_NP and cbut_bont,
> because queries relaxed no further than the GENUS. Relaxing to the FAMILY
> (Filoviridae) is what makes a profile buildable.

That machinery was written, committed, and never executed. `fetch_protein()`
walked `aa_queries()` and `return`ed at the first rung yielding **any** record:

```python
for query in aa_queries(row):
    uids = ez.esearch("protein", query)
    if not uids: continue
    ...
    return records[:want], query     # <- first non-empty, not first sufficient
```

Rung 1 is `NP[Gene Name] AND "Zaire ebolavirus"[Organism] AND
srcdb_refseq[Properties]`. RefSeq carries exactly one protein per locus per
virus species, so rung 1 returned 1 record, always, and the family rungs at the
bottom of the ladder were never issued. The fix for the starvation was sitting
below the line that made the starvation inevitable. `--dry-run` printed only
`aq[0]`, so no plan ever showed the rungs that did not run.

**The obvious fix is a false-positive generator, and it was measured, not
reasoned about.** Accumulating across the whole ladder until `want` records are
in hand does make every HMM buildable — and it populates markers with the wrong
organism, because the ladder's tail leaves the organism on purpose:

| marker | class | what pooling the whole ladder produced |
|---|---|---|
| `henipa_cedar` | `absence_expected` | **8 of 8 Measles virus** nucleoproteins, 0 Cedar |
| `nipah_L` | `required` | **7 of 8 non-Nipah**: Langya, Achimota pararubulavirus 3, *Respirovirus rupicaprae*, 2x Parajeilongvirus, a bat morbillivirus |
| `rvfv_Gn` | `required` | 8, of which 1 is RVFV; the rest came from the clade rungs |
| `cbut_bont` | `discriminating` | 2 real BoNT/E plus a thiazole biosynthesis protein, an MFS transporter and a 3-hydroxybutyryl-CoA dehydrogenase |

The last one is not a clade problem. `bont[All Fields] AND "Clostridium
butyricum"[Organism]` matches the *strain designation* of "C. butyricum E4 str.
BoNT E BL5262", so every protein that strain deposits is a hit.

This matters because `marker_evidence.py` states the scoring rule plainly: *"A
marker is PRESENT if ANY applicable tier passes."* These are not diluted
profiles. A measles-positive sample would have made `henipa_cedar` PRESENT, and
`henipa_cedar` is the `veto_if` for **both Nipah and Hendra** — the finding
would have been stood down by a childhood vaccine-preventable virus. `nipah_L`
is `required` for Nipah; a bat paramyxovirus would have satisfied it.

**So the ladder is split by what a rung can be trusted for.** `aa_query_plan()`
returns `(query, poolable)`, and a rung is poolable only if it stays inside the
declared organism **and** matches on a name field. Pooling runs over those; the
loose and clade-relaxed rungs survive only as a single-rung last resort, which
is exactly the pre-pooling behaviour, so a marker that only they can serve is
unchanged.

**And clade relaxation is now keyed on `marker_class`, not on
`query_strategy`.** It is allowed only for `pathotype`, where the gene is the
finding and the host species is not — BoNT in *C. butyricum*, Stx in *Shigella*
— which is the case `marker_rules.yaml` is built around and that D29 depends on
(`veto_if: cbut_bont` exists to ATTRIBUTE a cross-species hit, not suppress it).
For every other class the organism is the claim.

**This gate costs nothing that was working, and that is measured rather than
assumed.** Every `gene_family` row that was already `built` filled from
in-organism rungs alone, `cbot_bont` included (8 of 8 *C. botulinum*). An audit
of all 126 protein tiers found no pre-existing off-organism record that was not
a nomenclature synonym (*Orthomarburgvirus marburgense* for Marburg,
*Mammarenavirus choriomeningitidis* for LCMV, Poliovirus 1/3 under Enterovirus
C) or a legitimate multispecies RefSeq entry (Stx1 under Enterobacteriaceae).

**Result on the shipped pack.**

| | before | after |
|---|---|---|
| `built` | 132 | **143** |
| `partial` | 13 | **4** |
| `failed` | 5 | **3** |
| `curation_required` | 4 | 4 |
| `motif_todo` | 9 | 9 |
| HMM profiles in `markers.hmm` | 25 | **34** |
| rows in `curation_required.tsv` | 31 | **20** |

Nine HMM tiers now build from 8 in-organism sequences each: `ebola_NP`,
`ebolaS_NP`, `marburg_NP`, `marburg_GP`, `junin_GPC`, `rvfv_L`, `nipah_L`,
`nipah_G`, `eeev_nsP`. Two do not, and **the gate is why that is the right
answer**: `cbut_bont` has 2 genuine BoNT/E proteins and `rvfv_Gn` has 1, which
is precisely what D27 recorded ("only 2 protein sequences") before pooling
briefly papered over it with junk.

**Three retrieval failures were catalogue defects, not absent sequence.**

* `henipa_cedar` built **no query at all** and was recorded `failed`.
  **CORRECTION (D31): the fix described below is vacuous and this marker is
  still `failed`.** The same commit added the `nucleoprotein` alias, so the row
  now builds queries and never reaches the new branch; the static suite check
  then guarantees no shipped row can. The diagnosis in the rest of this bullet
  is sound, the remedy as shipped is not. Its gene
  symbol is `N`, one character, below `MIN_QUERY_NAME_LEN`, with aliases `-`, so
  `query_names()` returned `[]` and not one request was ever issued. `failed` is
  documented as "retrieval was attempted and returned nothing usable" — evidence
  that NCBI has nothing. This was silence from an instrument never switched on,
  filed as a negative result: the same conflation D26 removed from tier
  reporting. An unsearchable row is now `curation_required` with a loud
  per-marker error, and `tests/test_pipeline.py` has a static check that every
  retrievable row builds at least one query (mutation-verified: reverting the
  alias turns it red). Eighteen sibling rows carry a one-letter gene symbol
  rescued by an alias; this is what stops a nineteenth arriving without one.
* `bthai_ara` searched `araBAD[Gene Name]`, an **operon** name NCBI does not
  index: 0 hits, while `araA`/`araB`/`araD` return 2/2/15. Now `built`. It is
  the `veto_if` for **both** *B. pseudomallei* and *B. mallei*.
* `vchol_wbeO1` searched the **cluster** name: `wbeO1` and `rfbO1` both 0, while
  `wbeT`/`rfbV`/`wbeW` return 40 each. Now `built`.

**"Found nothing" and "found forty genomes and discarded every one" used to
produce the same note.** `fetch_nucleotide` now returns a third element
explaining an empty result, and says which filter did it. `mpox_C3L` reads
*"3 queries tried; they MATCHED 40 record(s), but none was usable as a locus"* —
195 kb MPXV genomes whose CDSs are named `MPXV-UK_P2-0NN`, not `C3L`. That is a
naming mismatch needing an ORF designation, and it had been pointing curators at
the query.

**Four gaps remain, each now carrying its measured diagnosis in the catalogue
note so it is not re-derived:**

| marker | why it is still open |
|---|---|
| `henipa_cedar` | Cedar's N is deposited as bare `N` under *Henipavirus cedarense* (3 nuccore hits); nothing is annotated "nucleoprotein" under either name. Needs an accession, not a query. |
| `bruc_S19_eri` | **Must not be "fixed" by retrieval.** It is a DELETION. `eryC` (19) and `eryD` (3) retrieve easily, but their presence is the wild-type state, and feeding them here would make `vaccine_if` fire on wild-type *B. abortus* and stand down a real select-agent detection. Needs absence-region logic, which does not exist. |
| `orthopox_ATI`, `mpox_C3L` | ORF/OPG designations needed as aliases. Guessing one would put the wrong gene into Variola's `confirm_requires`. |
| `lmono_LIPI4` | LIPI-4 is CC4-specific and absent from EGD-e, so the `lmo*` tags are the wrong reference. Needs a curated CC4 accession. |

Eleven new self-test checks on the ladder split and the relaxation gate, plus
one static suite check. `--dry-run` now prints every rung and names the
3-sequence minimum, because a plan that showed only rung 1 is how a dead
relaxation stayed invisible for its whole life.

---

## D31 — The protein tier's coverage gate was unreachable for 107 of 125 markers, and a failed Entrez request was being written as a verdict about NCBI

2026-10-15. Two findings from a Fable review pass, both in the shape the
standing principle names: a result indistinguishable from a clean negative.
The first is partly a regression D30 introduced.

### The protein coverage denominator

`fasta_lengths()` returned `marker_id -> SUMMED length` over every allele a
marker holds, and `aa_evidence()` divided the union of target intervals by it.
A marker carrying 8 alternative representations of one locus therefore capped
at **0.125** however well any single allele was covered, against catalogue
floors of 0.50-0.60.

**This is the third place the same defect has been found.** `nt_evidence()`
fixed it and wrote down why:

> BREADTH IS PER SEQUENCE, NOT OVER THEIR SUMMED LENGTH. A marker often carries
> several ALTERNATIVE representations of one locus [...] Dividing covered bases
> by the summed 1878 bp caps breadth at 0.5 even when one allele is covered end
> to end, so a min_coverage of 0.80 is structurally unreachable at ANY depth.

It was never carried to the protein tier. Measured on the shipped pack: **107 of
125** protein markers could not reach their own `min_coverage` on a perfect,
full-length, high-identity hit. `ebola_NP`, `nipah_L`, `marburg_NP`, `lassa_NP`,
`junin_GPC` and most of the other divergent-agent markers were among them — the
tier that exists to see through synonymous divergence was dead for exactly the
markers it was built for.

The harm is not a missing call. The row writes a confident **ABSENT**, which
feeds `absent()` as TRUE and can fire a `veto_if` or `vaccine_if` rule.

**D30 made nine of them worse.** Raising `ebola_NP` and eight siblings from 1-2
alleles to 8 cut their best achievable coverage from ~1.0 to 0.125. D30 gained
nine HMM tiers and killed nine protein tiers doing it, and said nothing about it
because nothing measured it.

**The fix mirrors `nt_evidence()` exactly:** lengths are keyed per reference
sequence, coverage is computed per target and the best one decides. Identity,
bits and e-value now come from the best-scoring hit **on that same sequence**,
so the two numbers the gate compares describe one protein — they previously came
from the best hit anywhere in the marker while coverage was unioned over every
allele, so a marker could pass identity on one allele using coverage accumulated
from another.

Measured after, against the real `markers.aa.faa` with a synthetic perfect
full-length hit on each marker's longest allele: **125 of 125 reach coverage
1.0000 and state PRESENT**, where 107 previously could not pass at all. The gate
still rejects — the same hits truncated to 20% give **125 of 125 ABSENT**, and
30% on each of four alleles does not sum to a pass.

Two new columns, `aa_n_seqs` and `aa_best_ref`, record which allele the
denominator came from. `aa_n_seqs` counts the SEARCH INDEX, not the manifest:
`concat_files()` collapses near-identical alleles when it builds
`markers.aa.faa`, so the manifest's `n_aa_seqs` (the archive) and this column
disagree by design and the gap is large — ebola_NP is 8 archived, 2 searched.
The coverage denominator needs the searched number.

The protein and HMM tiers also had **no reference-naming guard**, which
`nt_evidence()` has: a target present in the hits but absent from the FASTA gets
no length, so its coverage computes as 0 and the marker reads ABSENT in silence.
Now warned.

Eight `aacov` checks, mutation-verified: reverting the denominator alone turns
four of them red. The suite passed throughout the defect's life and had nothing
pinning this.

### A failed request is not an empty result

`Entrez.esearch()` and `efetch()` catch `EntrezError`, log at `LOG.debug` and
return `[]`/`""`, so the two `except EntrezError` handlers in `process_marker()`
are **unreachable** and a rate-limit block, a 5xx exhaustion or a truncated
efetch was recorded as `failed` — documented as *"retrieval was attempted and
returned nothing usable"*, i.e. a verdict about what NCBI holds. D30's own
`_tried_phrase` then asserted *"N queries tried, all empty"*, a positive claim
about a search that did not complete, and the marker landed on the curation list
sending a human to look for sequence that may well exist.

This host was HTTP 429'd throughout the 2026-10-15 rebuild, so it is not
hypothetical.

**A sixth status, `retrieval_error`**, now separates "the search did not
complete" from "the search found nothing". Failures are counted on the client
and logged at WARNING, `process_marker()` compares the count before and after
each marker, and a `partial` row whose retrieval also errored says so rather
than presenting its gaps as established. The run **exits non-zero** when any
request failed, because a build throttled end to end otherwise exits 0 and
leaves a manifest that reads as a measurement. `retrieval_error` also joins
`failed` in the blocking check for required high-consequence markers — shipping
a pack built from one is worse than shipping a known gap, because the gap is
real and the error is an artefact of the network.

203 passed, 0 failed, 2 skipped.

**Still open from the same review pass, not addressed here:** `flua_HA_H5` and
`flua_HA_H7` hold identical pan-subtype reference sets (1 of 7 is H5, 1 of 7 is
H7); `C3L` matches `SUBSTITUTION_RE` and is dropped before any query is built,
so `mpox_C3L` was filled from the `[All Fields]` fallback with B5R and
haemagglutinin; `nipah_N` and `rvfv_N` carry L polymerase and phosphoprotein;
132 of 163 manifest rows predate D30 and the resume drift-check compares neither
`gene` nor `aliases`, so no rerun refreshes them without `--force`; D30's
`curation_required` branch for unsearchable rows is dead code, because the same
commit added the alias that stops `henipa_cedar` reaching it — the commit
message and D30 both claim otherwise and are **wrong**.

---

## D32 — The organism was right and the subtype was the claim, and a poxvirus ORF name has the shape of an amino-acid substitution

2026-10-15. The remaining wrong-sequence findings from the Fable pass, and the
rebuild that applies them. Both are the same kind of error as D30's clade
relaxation — a reference set that cannot support the claim the marker makes —
reached by two routes the organism gate does not see.

### A third kind of scope

`clade_relaxation_allowed()` asks one question: does this rung stay inside the
declared organism. For `flua_HA_H5` the organism is right and the **subtype** is
the claim, so that question waves through a reference set exactly as wrong as
the clade rungs it blocks.

Measured on the shipped pack: `flua_HA_H5` and `flua_HA_H7` held the **same
seven HA proteins** — H5N1, H1N1, H3N2, H1N1, H2N2, H9N2, H7N9. One of seven is
H5; one is H7. Both rows were `built` with `has_hmm=yes`, so a profile was built
over the mixture. They are the two markers for highly pathogenic avian
influenza, and `marker_rules.yaml` lists them under `report_always`, so any
seasonal influenza A sample put "H5 marker present" and "H7 marker present" in
front of an analyst.

The winning rung for both was `HA[Gene Name] AND "Influenza A virus"[Organism]
AND srcdb_refseq[Properties]`. The `H5`/`H7` aliases contributed nothing,
because the gene symbol's rungs come first and already satisfied `want`.

**The fix is a declared scope, not an inferred one.** A new
`query_strategy=subtype_in_organism` says the aliases are tokens a retrieved
record must carry. NCBI encodes the subtype in the organism string — *"Influenza
A virus (A/goose/Guangdong/1/1996(H5N1))"* — so it is checkable **on the record**
rather than guessable from the query, which is the distinction the standing
principle's first lesson draws. The filter runs inside the pooling loop, so a
rung that returns nine wrong subtypes and one right one contributes the one and
the search keeps going.

`H5` matches `H5` and `H5N8`; it does not match `H15N1` or `H51`. The strategy
also removes these rows from the `gene_family` branch, so they no longer get
clade relaxation either — being `pathotype`, they were entitled to it.

### A poxvirus ORF name is shaped like a substitution

`SUBSTITUTION_RE` (`^(TR\d{1,4}|[A-Z]\d{1,4}[A-Z*]|\d{1,4}[A-Z])$`) exists for
alias lists like `cyp51A`'s `TR34;L98H;TR46;Y121F;T289A` — resistance
substitutions that are not searchable gene symbols. **Every value it was written
for is an alias.** Applied to the `gene` column it ate real symbols:

| row | gene | why it matched |
|---|---|---|
| `mpox_C3L` | `C3L` | poxvirus ORFs are `<fragment><n><L\|R>`, i.e. `[A-Z]\d{1,4}[A-Z]` |
| `pfal_18S`, `ehist_ssu`, `crypto_ssu`, `nfow_18S`, `bmand_18S` | `18S` | matches `\d{1,4}[A-Z]` |

`mpox_C3L` is the one that was harmed: its declared symbol was dropped before
any query existed, the ladder fell through to the `[All Fields]` last resort,
and the marker was built from **B5R, an EEV membrane glycoprotein and
haemagglutinin** — not the complement control protein, and useless for the clade
I/II call it exists to make. Status `partial`, i.e. testable.

The five `18S` rows are all `required`, and all survived on their `SSU rRNA`
alias and the rRNA title fallback — degraded rather than broken, searched
without their most precise term. The reviewer found `C3L`; the `18S` rows came
out of checking the regex against the whole catalogue rather than the one row.

**A curator who writes a symbol in the `gene` column means it is the gene**, so
the filter now applies to aliases only. The length floor still applies to both:
a one-character name is unsearchable whatever column it came from.

### The rebuild

The 2026-10-15 manifest carried two populations — 31 rows rebuilt under D30's
gate and **132 still produced by the pre-D30 code**, because the resume
drift-check compares `taxon`, `marker_class`, `evidence_type` and
`query_strategy` and none of those changed. It compares neither `gene` nor
`aliases`, which are exactly the fields D30's and this entry's remedies alter,
so `--force` is the only way to refresh them. Those stale rows still held what
the gate exists to stop: `nipah_N`, a **`required`** marker for Nipah, carried
an RNA-directed RNA polymerase L and two phosphoproteins among its seven
sequences, and `rvfv_N` carried the L polymerase inside the N marker. Both HMMs
were built over the mixture.

So this is a full `--force` rebuild of all 163 markers, which is also the first
one where a throttled run cannot be mistaken for a result (D31).

### The rebuild, measured

All 163 rows now carry `retrieved_utc 2026-10-15`, so the two-population
manifest is gone. Status totals are unchanged — 143 built, 4 partial, 3 failed,
4 curation_required, 9 motif_todo — and the three remedies landed:

| marker | before | after |
|---|---|---|
| `mpox_C3L` | `partial`, nt=0, built from B5R and haemagglutinin | **`built`**, a real `[gene=C3L]` CDS plus two complement control proteins |
| `flua_HA_H5` | 7 pan-subtype HA | **3, all H5N1** |
| `flua_HA_H7` | the same 7 | **3, all H7N9/H7N1** |
| `nipah_N` | 7, incl. an L polymerase and two phosphoproteins | **4, all Nipah nucleoprotein** |

`hendra_N` went `built` -> `partial` (3 -> 1 protein, HMM lost). That is an
honest loss, not a regression: it is `required` class, so D30's gate denies it
clade relaxation, and its previous three included off-organism sequence.

**A `--force` rebuild is not reproducible, and that is now recorded rather than
discovered twice.** Fifty rows changed sequence counts, in both directions —
`cocci_Ag2` 1->8, `toxo_GRA6` 1->8, `sars2_N_E` 2->8, but also
`ypest_ppcp1_pla` 2->1, `bpert_IS481` 3->1, `cchf_L` 8->3. `--force` bypasses
the Entrez cache, so every rung is re-ranked by NCBI as it stands today, and the
first productive rung's contents differ from the one cached in September. This
project pins its taxdump for exactly this reason; the marker database has no
equivalent pin, so two builds a month apart are not comparable artifacts.

**The 14 lost requests cost nothing, and establishing that took a hand
correlation that should not have been necessary.** They fell on `bthai_ara` (8),
`shig_stx1_Sd1` (2), `pjiro_mtLSU`, `crypto_ssu`, `bmand_mito` and `bmand_18S`;
all six recovered from other rungs and came out `built` with full sets, and none
of them is among the fifty rows whose counts moved. So the thinning above is
NCBI variability, not throttle damage — but the only way to know that was to
`awk` the log against the manifest.

A `built` row is trusted and SKIPPED by the resume path, so a row built while
the service was refusing requests cannot be told from a clean one. A new
`retrieval_errors` manifest column records the per-marker count, the run end
names every affected row — `NOT ESTABLISHED` for a `retrieval_error`, `SUSPECT`
for one that recovered — and prints the `--only-marker` arguments to refresh
them. Without it, "re-run the ones that were affected" is not an executable
instruction, which is the recovery plan `retrieval_error` assumes. The column is
additive: a manifest predating it reads as blank rather than as zero errors, and
`csv.DictWriter` carries old rows forward unchanged.

The exit code stays deliberately blunt — any failed request fails the run, even
when every marker recovered, as happened here. The alternative considered was to
fire only when a marker ended non-`built`, and the data rejects it: a marker can
end `built` and short, which is the case that most needs flagging and the one
that refinement would have hidden.

## D33 — An arm that was enabled and died looked exactly like an arm that was switched off, and the one signal that survives a mislabelled novel agent never reached the page

**Found 2026-10-15, by reading a CAVEAT the code had written about itself.**
`modules/local/detect.nf` carried this above METABULI's `errorStrategy 'ignore'`:

> screen_summary.tsv is the record, and it reaches neither SCORE_EVIDENCE nor
> MAKE_REPORT — it is only published. [...] if this arm were the ONLY one that
> would have raised a taxon, its ignored failure loses that taxon with no trace
> in the report at all. The summary is published beside the results as the
> compensating record.

The compensating record did not compensate, for two independent reasons, and the
second is worse than the first.

### 1. The record could not record the thing it existed for

`fillMissing()` pads an absent arm with `[]`. That was the right fix for the
basename collision it replaced (two padded arms sharing `assets/NO_FILE` killed
PARSE_SCREEN with "input file name collision"), and it has a side effect nobody
wrote down: by the time `parse_screen.py` runs, **three different events are one
state**. A disabled arm, an arm that crashed under `errorStrategy 'ignore'`, and
an arm that emitted an empty report all arrive as the same missing `--<arm>-report`
argument.

The summary loop then made it worse by skipping what it could not measure:

```python
for label, table in (("krakenuniq", ku), ("kraken2", k2), ("metabuli", mb)):
    if not table:
        continue
```

So the failure the CAVEAT describes produced **no row at all** — and neither did
a disabled arm. The file whose entire job was to be the trace of an ignored
failure recorded that failure as absence, which is the same thing it recorded
for a deliberate configuration choice. `if args.out_summary and summary:` closed
the circle: on the one run where every arm died, the file was not written at
all.

Only the workflow knows which arms it asked for, so it now says so
(`--expected-arms`), and the four states are named: `ok`, `empty_report`,
`no_output`, `not_run`. Every arm gets a row whether or not it misbehaved,
because an answer printed only when something is wrong cannot be told apart from
a missing table.

### 2. The unclassified fraction was computed, logged, written, and never read

This is the larger loss, and it was not in the CAVEAT. `unclassified_fraction()`
carries a fourteen-line docstring explaining that it exists because of Perseus
(Bioinformatics 2026, doi:10.1093/bioinformatics/btag687): sequence from a family
absent from the database is not left unclassified, it is **confidently assigned
to the nearest known taxon**, approaching a 50% false-assignment rate on long
reads. The docstring's own words are that the unclassified fraction "is the only
signal that survives that failure mode, so it is reported on every sample rather
than thrown away."

It was thrown away. It went to `screen_summary.tsv` and to stderr, and
`screen_summary.tsv` reached nothing. An analyst reading the HTML report — the
artifact this project treats as the deliverable, on the stated principle that "a
caveat belongs on the finding and in the HTML, not in terminal output nobody
keeps" — had no way to see it. The shape is the one `roadmap.md` §0b calls THE
novel-agent failure mode, and the mitigation was built, tested, and then not
plugged in.

`--high_unclassified` compounded it: `parse_screen.py` has had the argument since
the function was written, and **the process never passed it**, so the 0.40
threshold could not be changed without editing the script. It is now a parameter,
and documented as uncalibrated — the right value depends on matrix and database,
and a clinical sample against PlusPF is not comparable to an environmental one.

### Why this does not reach SCORE_EVIDENCE

It is coverage, not evidence, and it must not move a call. A high unclassified
fraction is a reason to look harder at a sample, not a reason to raise or lower
a tier on any taxon in it; a dead arm is a hole in what was searched, not an
argument about what was found. Routing it into the scorer would make "the
metabuli job crashed" into a term in a tiering decision, which is the wrong kind
of fact in the wrong place. It goes to MAKE_REPORT and stops there.

### Addendum, same day: `not_applicable` is not a shade of `not_run`

Caught by the restarted `zymo_ont` run, on the first real ONT data this coverage
table has ever seen, and not by anything in the suite.

Bracken does not run on long reads — its model takes one fixed read length —
so the workflow drops it from the expected set for a non-Illumina sample. Dropped,
it fell through to `not_run`, whose note reads "disabled for this run", **on a
sample whose kraken2 row read `ok` and where nobody had disabled anything.** A
reader would reasonably conclude an operator had passed `--skip_kraken2`.

This is a smaller version of the same mistake the entry above is about: two
different answers to "why was this not searched" collapsed into one, where the
collapsed answer is the misleading one. An operator's choice and a property of
the platform are not interchangeable. `--inapplicable-arms` now carries the
second, and naming an arm as both expected and inapplicable is a contradiction in
the workflow's own description of the run and is rejected rather than resolved to
either.

Worth noting about how it was found: the suite passed, the unit tests passed, and
the preview passed. It took real data on the one platform the pipeline had never
processed. That is the same lesson `benchmarking.md` already records, arriving
again.

### Bracken is deliberately exempt from the alarm

It is the one arm whose absence is not worth flagging, and the reason is not
that it matters less. It never appears in `support`, so it cannot raise a
candidate and no taxon can be lost by its dying; and it **exits non-zero by
design** when a report has no reads at the requested level, which is the correct
outcome for a clean negative control. Flagging that would put a warning on every
spotless control in the run and teach the reader to skip the column. It still
gets a row — the question the table answers is what was searched — with a note
saying the absence is routine.

### What the report says now

A dead enabled arm and a `HIGH_UNCLASSIFIED` sample each raise a banner above
the findings, where a reader meets them before deciding the page is clean, and
both are detailed in a **Stage-1 screen coverage** table that lists every arm for
every sample with the failures sorted to the top.

### Verified

`tests/test_pipeline.py` gained a `report` group — **make_report.py had no test
coverage at all before this**, which is part of why a whole input could be
missing from it unnoticed — and eight new assertions in `screen`. They pin the
distinctions that are easy to lose again: an arm nobody asked for is `not_run`
and never a failure; an enabled arm that produced nothing is `no_output` with
`ARM_PRODUCED_NOTHING`; bracken's absence is recorded without a flag; an arm with
no read ledger reports **blank, not `0`** (`0` unclassified of `0` total is a
claim about the data, and sylph has no read ledger while bracken's would be
kraken2's counted twice); a run where every arm died still writes its table; a
misspelled `--expected-arms` value fails rather than becoming a silently
unenforced expectation; and a report rendered without a summary omits the section
rather than drawing an empty one that reads as all-clear.

Beyond the suite: `nextflow run -preview` on the plumbing samplesheet to check
the channel arity, and two sylph-only end-to-end runs on the plumbing positive
control — deliberately off the 123 GB PlusPF index, which would have put the box
back under the `systemd-oomd` pressure trap recorded in STATUS.md.

### The failure, reproduced end to end rather than argued

The second run pointed `--sylph_db` at a file sylph cannot read, so the only
enabled arm died **through the real channel topology**: `SYLPH_PROFILE` exited 1,
`errorStrategy 'ignore'` swallowed it, `fillMissing()` padded it to `[]`, and the
run continued to MAKE_REPORT and exited having analysed nothing.

The input is the plumbing positive control: 2M reads carrying *Y. pestis* at
5,000 RPM, whose scenario file says any outcome other than a confirmation is a
configuration fault.

| | old wiring | now |
|---|---|---|
| `findings.tsv` | header only, 0 rows | header only, 0 rows |
| report | every tier renders "None."; the only banner is the unrelated "No negative control" | banner: *"A stage-1 classifier was enabled and produced nothing: sylph on 1 sample(s)"* |
| mention of the dead arm anywhere in the HTML | **none** | banner + a `Stage-1 screen coverage` row, `no_output` / `ARM_PRODUCED_NOTHING` |
| `screen_summary.tsv` | not written — `if summary` skipped the arm that produced nothing, leaving nothing to write | 5 rows, 4 `not_run` and 1 `no_output` |

Re-rendering that same `findings.tsv` with the summary withheld reproduces the
old page exactly: zero occurrences of "sylph", of "coverage", or of any word
about an arm having failed. A clean screening report on a sample with a
notifiable organism in it at 5,000 RPM, and the run exited reporting success.

This is the same shape as the stage-2 read-count filter that silently discarded
every sylph hit, and as the thirty-odd before it: *produces an empty or degraded
result indistinguishable from a clean negative*. It is also the one case the
CAVEAT in `detect.nf` predicted in writing and then mitigated with a file that
nothing read.

## D34 — Twelve references had never been masked, and the record of what was examined was keyed on a taxid that owns twelve of them

**Found 2026-10-15 while implementing a different fix**, which is the only reason
it was found at all. STATUS.md carried a to-do: the `SHARED_SEQUENCE` flag is
derived from run-time contig collisions, so a `--only-taxid` run of one partner
produces correct numbers with an empty flag, and "the durable fix is to derive
sharing from the refs manifest's `accession` column, which is an artifact rather
than run state."

Doing that required looking at the `accession` column, and the column says
something the to-do did not anticipate.

### The human pack has eight shared assemblies, not two

| accession | entries |
|---|---|
| GCF_000833355.1 | 263 *F. tularensis* / 264 *F. novicida* |
| GCF_022354085.1 | 622 *S. dysenteriae* / 620 *Shigella* |
| GCF_000848505.1 | **11266 Filoviridae** / 186538 Zaire ebolavirus |
| GCF_000854085.1 | **11266 Filoviridae** / 186539 Reston ebolavirus |
| GCF_000855585.1 | **11266 Filoviridae** / 186540 Sudan ebolavirus |
| GCF_000857325.3 | **11266 Filoviridae** / 3052505 Marburg marburgvirus |
| GCF_000889155.1 | **11266 Filoviridae** / 565995 Bundibugyo virus |
| GCF_000943645.1 | **11266 Filoviridae** / 3052505 Marburg marburgvirus |

Six of the eight are the Filoviridae group entry (D24/D25) sharing an assembly
with the member taxon it was built from. The run-time detection could never have
seen any of them: the partners have **different taxids**, so no `--only-taxid`
invocation puts them in one run, and a full-pack run is the only arrangement in
which they collide at all.

### And then: twelve references had no masking result whatsoever

Reconciling the refs manifest against the masked directory by hand — 371
references, 359 BEDs — all twelve of `11266_Filoviridae*` were missing. The
group entry was added to the pack on 2026-10-13/14; the masking pass ran
2026-10-03/04 and nothing re-ran it.

Six of the twelve share an assembly with a member taxon, so their contigs appear
in that member's BED and they were masked by accident. **The other six carry
contigs that appear in no BED anywhere in the pack:**

| accession | contig | in any BED |
|---|---|---|
| GCF_000888475.1 | NC_014372.1 | no |
| GCF_000896415.1 | NC_016144.1 | no |
| GCF_003505815.1 | NC_039345.1 | no |
| GCF_004788975.1 | NC_055175.1 | no |
| GCF_004788995.1 | NC_055176.1 | no |
| GCF_013088285.1 | NC_055510.1 | no |

`coverage_stats.load_mask()` reads every BED in the directory into one
contig→intervals map. A contig that is in no BED is simply absent from that map,
which is exactly what a contig with no masked intervals looks like. So those six
references were scored with their full length in the breadth denominator and no
`reads_in_masked_regions`, and **nothing distinguishes that from "examined and
found clean."** This is D5's lesson one level up: not an unmasked interval but an
unmasked reference.

The practical cost here is probably small — every filovirus reference that WAS
masked came back at 0 bp, and these are 19 kb RNA virus genomes — but "probably
small" is a statement that could not be made from the artifacts, only from
re-running them. The whole point of the mask manifest is to not have to.

It also lands on the one entry that exists to catch a divergent agent. D25
measured the Filoviridae group entry recovering Bundibugyo 3/3 at zero
specificity cost; its references were the unmasked ones.

### The record was keyed on the wrong thing, in two places

`mask_manifest.tsv` is merged, not rewritten — a fix already made here once,
after an interrupted run destroyed the record of 296 masked genomes. The merge
key was **taxid**:

```python
fresh = {str(row["taxid"]) for row in summary}
... if old["taxid"] not in fresh: carried.append(old)
```

A taxid is not a reference. 11266 owns twelve references and 186538 owns three,
so masking one of 186538's three files writes one fresh row and drops the other
two: a taxid with three BEDs and one manifest row, exit code 0. The same
truncation the comment above it says was fixed, surviving one level down because
the key was never the thing being written. It is now `bed_path`, which is unique
per reference and was already a column.

`tests/test_pipeline.py` had the identical bug, and it is worth recording because
it is why the suite was green. The `shared_assembly` group read the manifest into
`{r["taxid"]: r}` — so eleven of 11266's twelve rows were discarded by the
comprehension. The pairs it checks are one-reference taxa, so it passed. It was
asserting less than it appeared to.

### What is now in place

`ref_stem(row)` is the single definition of a reference's identity, and the
docstring says what it is not. `SHARED_SEQUENCE` is derived from `accession`, so
it no longer depends on how the script was called. And `reconcile()` compares
**every** row of the refs manifest against the results on disk, on every exit
path including `--dry-run` and including the "nothing to do" early return — which
is the run most likely to be read as "the pack is complete", and the one where a
`--only-taxid` reference never enters `todo` at all. It names each reference with
no result, prints the command that would fix it, and returns 1. Blunt on purpose,
the same way D32's retrieval exit code is: the run did what it was asked and the
PACK is still incomplete, and an exit code reflecting only the former is how this
sat unnoticed for twelve days.

A `--self-test` was added (the script had none) covering the three conventions
offline, with no pack, genome or minimap2: that twelve references under one taxid
stay distinct, that a shared assembly is found from the manifest without reading
a genome, and that an **empty** BED counts as examined rather than missing —
because low-complexity masking finds something in essentially every genome, a BED
exists either way, and a later change that "fixed" emptiness into absence would
re-create the confusion this entry is about.

### Self-review of this change found two more, both of my own making

Worth recording because both are the same category as everything above.

**`reconcile()` compared full paths.** `bed_path` is written absolute, so a pack
copied to another machine, mounted elsewhere, or masked through a symlink would
have matched nothing and the function would have reported all 371 references
unmasked. A false alarm in a check whose only value is being believed is worse
than no check: the habit it teaches is to stop reading the warning. Compared by
basename now, and the self-test writes a manifest under a foreign root and
requires it to reconcile.

**MAKE_REPORT could be killed by its new input.** It sets no `errorStrategy`, so
an exception there ends the run **after** stage 1, stage 2, the marker tier and
scoring have all finished — trading the result for the view of it. That is
precisely what bf9744d fixed in PAVIAN_EXPORT one commit earlier, and wiring a
new file into the last step of the run re-opened it: `read_tsv` raises on a
truncated or non-UTF-8 summary, and `float()` raises on a non-numeric
`unclassified_fraction`. It now degrades to a note **that says the coverage is
unknown** — not to silence, because an omitted section is what a run with no
summary at all produces, and a reader cannot tell "no arms to report" from
"could not read what the arms did".

### The run says it too, not only the builder

`mask_pack_refs.py` reconciling itself only helps someone who re-runs the masker.
A pipeline run uses whatever is on disk, and the workflow's own check was
all-or-nothing — the `masked/` directory exists or it does not — so it could not
see the state the pack was in for twelve days. `warnIfPackPartlyMasked()` now
compares the refs manifest to the BEDs present and names what has no result.

It WARNS rather than exits, unlike the GTDB-namespace guard beside it. A
partially masked pack still produces findings worth having; a GTDB-namespaced
index produces an arm that contributes literally nothing. Those deserve different
answers.

Verified both ways on the real pack — silent when complete, and `1 of 371 ...
(Shigella)` with one BED moved aside and restored byte-identical afterwards. The
first attempt at that verification printed nothing and looked like a pass; it was
a Groovy compile error (`considered++` is rejected by Nextflow's strict script
syntax), which is worth recording as its own small lesson: a check that produces
no output is not the same as a check that passed.

### How it was found, which is the part worth not repeating

By diffing a directory listing against a manifest, by hand, while doing something
else. Nothing in the pipeline, the suite, the mask manifest or the run logs said
371 ≠ 359. D32 closes with almost the same sentence about a hand `awk` that
should not have been necessary; this is the second time, and the remedy is the
same both times — have the artifact reconcile itself and fail when it cannot.

## D35 — The first ONT measurement: specificity is clean, every call is capped by one provisional threshold, and the instrument scored a false FAIL on the only organism it got right

`runs/zymo_ont` completed 2026-10-16 after 15h 36m, 425 CPU-hours. It is the
first ONT data this project has ever processed end to end, and the only evidence
the D10 identity floors have ever had.

### It died at 4 h last time for exactly the reason D18 said

`MARKER_SEARCH_AA` on `zymo_log_ont_ERR3152366` took **7h 20m 28s** at 2262% CPU
(22.6 cores). The old 4 h base default was not marginally short; it was off by
nearly a factor of two, and the run would have died at the same place. Everything
else in the run is minutes — COVERAGE_STATS 2–19 m, MARKER_ALIGN_NT 3.5 m,
MARKER_SEARCH_HMM 12.6 m, PAVIAN_EXPORT 12 s. One process is 97% of the wall
clock, which is also why its `maxForks = 1` matters.

### Specificity: clean

84 taxa on the truth-absent list. Four produced findings, all held at
`attributed_to_neighbour`, none reaching `supported` or above:

| taxon | breadth | identity | attributed to |
|---|---|---|---|
| *Bacillus anthracis* | 0.0825 | 0.7010 | taxid 1405 — *B. subtilis* IS in the standard |
| *Shigella dysenteriae* | 0.7810 | 0.8203 | 624 — *E. coli* is present, and Shigella is a pathotype within it |
| *Candida auris* | 0.0053 | 0.7339 | 45357 |
| *Cryptococcus gattii* | 0.0012 | 0.7490 | 40410 — *C. deneoformans* |

**Zero false confirmations.** Tier PPV: `supported` 2/2 = 1.000.

### Sensitivity: one of four, and the cap is a single threshold

Four watchlist taxa are genuinely present. *L. monocytogenes* reached
`supported`; *E. coli*, *S. enterica* and *C. neoformans* were all absorbed by
`attributed_to_neighbour`, giving that tier a PPV of 0.556 (5 of 9 calls are on
taxa that are actually there).

The *Listeria* result is the one to read closely, because **the only gate it
failed was the identity floor**:

```
read threshold met; failed: identity 0.8579 >= 0.9000 (ont/default)
marker rule satisfied: lmono_hly and lmono_prfA
marker_verdict: confirm_requires = TRUE; support_requires = TRUE
```

breadth 0.9525 / 0.9646, rpm 136,951 / 821,645, neighbour_ratio 9.74 / 12.13,
both confirmatory markers PRESENT, `confirm_requires` satisfied. With any ONT
identity floor at or below 0.849 both calls are `confirmed_in_silico`. At 0.90
**no ONT sample in this run can reach `confirmed_in_silico` on anything** —
observed identities span 0.70–0.89 across all eleven findings.

### But this does NOT license lowering the floor, and the reason matters

The data is Nicholls et al. 2019 (GridION, ERR3152364/ERR3152366) — **R9.4.1**.
D10 reasoned the 0.90 default from **R10.4.1 simplex at ~97–98%**. Observed
0.85–0.89 is exactly what R9.4.1 should give. So the floor was not falsified;
**this benchmark cannot test it.** Fitting 0.90 down to R9 data would make the
gate far too permissive for the chemistry it is actually calibrated for, which
is the unsafe direction.

What the data DOES support, and independently of chemistry, is the pack's own
stated position that on ONT "read identity CANNOT do that job". Identity does not
separate true from false here — *E. coli* is present at 0.7999 while the absent
*S. dysenteriae* reads 0.8203. The gate is not discriminating; it is only
capping. Discrimination is being done by competitive mapping and the marker
tier, exactly as `pack.yaml` predicts.

**So the action is not a new number. It is that the ONT arm still has no data on
the chemistry it is calibrated for**, and `roadmap.md` F3 already records that no
independent R10.4.1 benchmark of ONT shotgun classifiers exists to borrow from.
Until R10.4.1 data is run, `confirmed_in_silico` on ONT should be understood as
unreachable by construction rather than as a high bar.

### The evaluator scored a false FAIL on the one organism the pipeline got right

Scoring this run against the full manifest gave 5 pass / 5 fail, including:

```
FAIL zymo_even/Listeria monocytogenes min_tier: expected >= supported, observed no finding
```

*Listeria* was `supported` in the ONT sample. The scenario `zymo_even` spans an
Illumina sample and an ONT one, `evaluate()` walks the MANIFEST rather than the
run, and the Illumina partner simply was not in this findings file — so it scored
`tier=None` → NEGATIVE, `_best()` returned "no finding", and the scenario failed.
Restricted to the samples the run covered: 6 pass / 4 fail, and those four are
substantive.

This is not an edge case. Bracken takes one read length per run, so this project
splits one mock community across `zymo_rl100`, `zymo_rl150` and `zymo_ont` **by
necessity**, and every such run scored against the whole manifest charges itself
for the samples it did not contain. The instrument used to judge everything else
was the thing producing the wrong answer, and it produced it in the worst
possible place.

`--samplesheet` now names the samples a run analysed; manifest rows outside it
are excluded from grading and counted, rather than scored as misses. Without it
the old behaviour is kept and warned about — narrowing to "samples that produced
findings" would silently reclassify a genuinely clean sample as "not run", which
is the same conflation in the other direction.

**Third occurrence of one lesson.** The samplesheet is authoritative for which
samples are CONTROLS, because a spotless control produces no rows. PARSE_SCREEN
must be told which stage-1 arms were enabled, because a crashed arm and a
disabled one both produce nothing (D33). And now: a sample that ran and found
nothing is indistinguishable, in a findings file, from one that never ran. Each
time, absence of output was read as evidence of absence, and each time only the
thing that launched the run could say which it was.

### Also measured

The unclassified fraction on real ONT against PlusPF is **0.000525** (EVEN) and
**0.000136** (LOG) — a floor, not a typical value: ten organisms PlusPF knows
well, with no host and no environmental diversity. It is the first time the
number has been measured on long reads here.

*C. neoformans* was not raised at stage 1 at all in the LOG community, where it
sits at the bottom of a log-distributed standard. That is a genuine miss, at an
abundance far below any LoD this pipeline claims.

## D36 — The first real Illumina re-score flagged four correct answers as failed arms, and the warning that fires on correct behaviour is the one the reader stops believing

**Status:** fixed. Found mid-run, in the half of `runs/zymo_rl100_rescore` that had
already completed, four days after D33 introduced it.

### What happened

D33 gave the stage-1 screen a coverage table, so that an arm which was *enabled and
died* could no longer be read as a clean negative. It distinguishes five states, and
two of them raise a flag that is bannered at the top of the report:

- `no_output` → `ARM_PRODUCED_NOTHING` — the arm was enabled and emitted no report
  at all. Under `errorStrategy 'ignore'` that is an invisible crash. This is the
  state D33 existed to catch and it remains exactly as it was.
- `empty_report` → `ARM_EMPTY_REPORT` — the arm emitted a report with no usable rows.

The second one was wrong for sylph, and the first real Illumina data to reach it said
so. Four of eleven single-organism libraries — *P. aeruginosa*, *S. aureus*,
*E. faecalis*, *L. fermentum* — carried `ARM_EMPTY_REPORT` on sylph. All four were
correct. Those four organisms are not on the watchlist and are not near anything on
it, and sylph profiles against `human_pack_refs.syldb`: **the pack, 371 watchlist
references and nothing else.** It emits one row per containment hit, so a
header-only profile is how sylph says *none of the watchlist is in this sample*. That
is the commonest correct answer the arm has.

So the report for this run would have opened with a banner reading "A stage-1
classifier was enabled and produced nothing: sylph on 4 sample(s) … its absence is
NOT a negative result", on four samples where the pipeline had answered perfectly.

### Why the arms are not alike

The distinction is not the parser and not the arm's importance. It is what the arm
searched:

| arm | database | parsed | zero usable rows means |
|---|---|---|---|
| kraken2, krakenuniq, metabuli | comprehensive | unfiltered | **broken report** — these emit rows for any input at all: root, unclassified, host |
| sylph | the pack (371 refs) | per containment hit | **the negative answer** |
| bracken | re-estimates over kraken2 | — | nothing lost; it cannot raise a candidate (D33) |

`raises` already existed as a per-arm property and already exempted bracken, but it
does not capture this: sylph *does* raise candidates, so a sylph crash must still be
flagged loudly. What differs is only whether an *empty* report is informative. That
is now a separate property, `empty_is_negative`.

The row is still written, still says `empty_report`, and now carries a note that
says why it is not an alarm — `ran against the pack and matched no watchlist
reference - this IS the negative answer, not a missing one`. Suppressing the flag
must not suppress the explanation: the row is what tells a reader the arm ran at all.

### Why this was worth stopping for

This project's standing concern is the bug that produces an empty result
indistinguishable from a clean negative. This is its mirror image, and it is not a
cosmetic one. `ARM_PRODUCED_NOTHING` and `ARM_EMPTY_REPORT` share a column, a banner
and a sort rank. A banner that fires on four correct answers in an eleven-sample run
is a banner an analyst learns to scroll past — and the next thing it says may be that
metabuli died and a divergent agent went unsearched. The same reasoning retired the
`reconcile()` false alarm in D34's self-review: *a false alarm in a check whose only
value is being believed is worse than no check.*

### The cost of the suppression, and where it is paid

Suppression opens a real hole. An empty `.syldb`, or one truncated at build time,
**loads without error** and returns zero rows for every sample. Per-sample that is
now indistinguishable from a cohort that genuinely carries nothing.

Across a run it is not indistinguishable at all, and the report is the only place
that holds every sample at once — `parse_screen.py` sees one sample and cannot know.
So the alarm moved rather than vanishing: `make_report.py` now banners an arm that
*searched more than one sample and matched nothing in any of them*, as a database
problem rather than a crash. An arm empty on some samples and `ok` on others is a
normal mixed cohort and stays silent — one `ok` row proves the database loaded.

This run is the case that distinguishes them: sylph `ok` on 7, `empty_report` on 4.
No banner, correctly.

A corrupt or missing `--sylph_db` was already covered and is untouched: it makes
sylph exit non-zero, so no file is produced, so the state is `no_output` and the
alarm is `ARM_PRODUCED_NOTHING`. That is how D33 was proved in the first place.

### Verified

Reproduced on the real profile, not a fixture —
`runs/zymo_rl100_rescore/sylph_profile/zymo_isolate_pseudomonas_aeruginosa_illumina_ERR2935853.sylph.tsv`,
a genuine header-only sylph output, run through `parse_screen.py` at `HEAD` and at
the fix:

```
BEFORE  flag column: [ARM_EMPTY_REPORT]   + a WARNING on stderr
AFTER   flag column: []                   + the explanatory note, no warning
```

Six new assertions, **234 passed / 0 failed / 2 skipped** (was 228):

- a pack-restricted arm's empty report is the negative answer, and is unflagged
- it still says so in the note rather than going quiet
- a *crashed* pack-restricted arm — file absent, not empty — is still
  `ARM_PRODUCED_NOTHING`, so the suppression cannot reach the crash
- an arm empty on every sample is bannered as a database problem
- and is not described as a crash, because it ran
- a mixed cohort raises no database alarm

### Note on the in-flight run

The flags are written by `PARSE_SCREEN`, which had already completed for all eleven
samples when this was found, so `runs/zymo_rl100_rescore`'s own report will still
carry the stale banner. Nothing in it is a call: the screen summary deliberately
reaches `MAKE_REPORT` and never `SCORE_EVIDENCE` (D33), so no tier can move. The
summaries and the HTML will be regenerated by direct script invocation after the run
completes, rather than by `-resume`, which would re-hash `bin/` and risk redoing
eleven hours of `MARKER_SEARCH_AA` to fix a banner.

## D37 — `confirmed_in_silico` can fire on real data, and the one remedy the tool prints for the gate still holding it back did not parse

**Status:** fixed. Found by scoring the five completed samples of
`runs/zymo_rl100_rescore` while the other six were still running.

### The measurement the re-score was launched for

No run in this project's history had ever produced a `confirmed_in_silico` on real
data — rl100 43 `supported`, rl150 6, rl150_rescore 1, zymo_ont 2, zero confirmed.
The motivating question was whether the top tier was reachable on real data **at
all**, and the pack recorded the suspected cause: *S. enterica* in its own
single-organism library at breadth 0.9137, occupancy 0.9378, neighbour_ratio 1928,
held at `supported` by `mean_alignment_length 99 >= 100` — a gate D20 replaced with
`min_alignment_fraction_of_read` and never re-ran against the data that motivated it.

**D20's fix holds.** The same library, same organism, scored with current code:

```
all confirmation gates passed: breadth 0.9137 >= 0.0100; occupancy 0.9378 >= 0.1500;
identity 0.9870 >= 0.9800 (illumina/default); depth_cv 0.42 <= 2.50;
mean_alignment_length 99 >= 61 (illumina 60% of 101 bp read); mapped_fraction 1.000 >= 0.50
```

`99 >= 61` where it was `99 >= 100`. Every alignment gate is now reachable on real
2x101 Illumina data.

**And it was still `supported`,** held by something that is not a threshold:

```
CAPPED AT SUPPORTED: no marker rule block exists for this taxon, so the veto,
vaccine and confirmation checks never ran. ... Add a block to marker_rules.yaml, or
declare `no_markers_available: true` there to state that none exist.
```

That cap is D17a and it is correct — alignment geometry cannot tell wild type from a
vaccine strain, nor a target from a near neighbour carrying a diagnostic deletion.
`marker_rules.yaml` has 58 blocks against 90 watchlist taxa, so for the other **32,
including *S. enterica*, `confirmed_in_silico` is unreachable by construction**,
whatever the data looks like. (D17a's docstring predicted exactly this organism.)

### The bug

The cap names a way out. Following it fails the run:

```
ERROR: cannot parse marker_rules.yaml: unknown rule key 'no_markers_available'
for 'Salmonella enterica'; allowed: ['confirm_requires', 'rationale',
'report_always', 'support_requires', 'vaccine_if', 'veto_if']
```

`score_evidence.py` reads `spec.get("no_markers_available")` and prints the remedy;
`marker_rules.py`'s `RULE_KEYS` allow-list never contained the key. The only
documented escape from the cap was unusable, so the 32 taxa were capped with no way
to lift it short of writing marker biology.

### Why nothing caught it

`tests/test_pipeline.py` has a case named for this exact behaviour, and it passed:

```python
rules = tmp / "rules_none.yaml"
rules.write_text("taxa:\n  Yersinia pestis:\n    no_markers_available: true\n")
...
got = score(tmp / "d17a", cov, markers=_empty_markers(tmp))   # <- `rules` never passed
record("d17a:the cap says why, and how to lift it",
       "no marker rule block" in got["reasons"]
       and "no_markers_available" in got["reasons"], ...)
```

The file is written and then never handed to `score()` — the helper hardcoded the
shipped pack — and the assertion is satisfied by the string appearing in the **cap
message itself**. It asserted that the cap recommends the flag, never that the flag
works. A test can name the behaviour it means to pin and pin nothing; `score()` now
takes a `rules=` override so it can actually be exercised.

### `bool("false")` is True, which is why the flag is parsed strictly

Every other value in `marker_rules.yaml` is a string. Had this one been stored the
same way, `no_markers_available: false` would have **granted exactly the permission
it denies**, and `flase` would have granted it too. For a flag whose only job is to
lift the cap that keeps an unmarked taxon out of the top tier, that is the unsafe
direction. It is now coerced to a real bool from `true/false/yes/no/on/off/1/0`, and
any other value is a hard parse error rather than a truthy string. Same reasoning as
the comment already in `score_evidence.py`: *silence must not grant permission; a key
in the pack must* — and neither must a typo.

### Lifting the cap must not remove the disclosure

Nothing was written when the declaration *did* apply, so a declared gap and a
fully marker-confirmed call would have been indistinguishable on the page. A
`confirmed_in_silico` reached this way now carries:

> MARKER TIER DECLARED UNAVAILABLE for this taxon (`no_markers_available: true`), so
> the veto, vaccine and marker confirmation checks did not run. This call rests on
> alignment geometry alone; a vaccine or attenuated strain, and a near neighbour
> carrying a diagnostic deletion, cannot be excluded by that evidence.

### Verified

With the fix, and the flag declared for *S. enterica* in a scratch copy of the rules:

```
BEFORE  attributed_to_neighbour 19   screened_only 9   supported 1
AFTER   attributed_to_neighbour 19   screened_only 9   confirmed_in_silico 1
```

**The first `confirmed_in_silico` on real sequencing data in this project.** It is a
true positive: *S. enterica*, notifiable, in its own single-organism library.

**That declaration is NOT committed, and must not be.** *Salmonella* has
well-established diagnostic targets (*invA*, *ttrRSBCA*), so
`no_markers_available: true` would be a false statement about the biology. What the
experiment establishes is that the mechanism works; the correct fix for this taxon is
a real rule block, which is curation work.

Suite: **239 passed / 0 failed / 2 skipped** (was 234). Six new assertions, including
that `no_markers_available: false` does not lift the cap and that a mistyped value is
refused rather than read as permission.

### Side measurement: the re-score collapses `supported`, in the right direction

Like-for-like on the five completed samples:

| | Oct 7 | 2026-10-19 |
|---|---|---|
| `supported` | 15 | **1** |
| `attributed_to_neighbour` | 7 | 19 |
| `screened_only` | 7 | 9 |

Fourteen `supported` calls became `attributed_to_neighbour`, the same direction
rl150 moved on its Oct 15 re-score (6 → 1). These are single-organism libraries, so a
*Salmonella* call in the *E. coli* library IS a cross-attribution artifact and
`attributed_to_neighbour` is the right answer: this is specificity improving, not
sensitivity lost. The two libraries containing nothing on the watchlist
(*P. aeruginosa*, *E. faecalis*) produced no call above `attributed_to_neighbour`.

Two genuine misses remain, both near-neighbour cases where the neighbour is arguably
the harder claim: *E. coli* in its own library is attributed to *Shigella*
(phylogenetically **inside** *E. coli*, and sharing a pack assembly per D34), and
*C. neoformans* is attributed to *C. gattii*, its sibling in the same species
complex. Full grading waits for the remaining six samples and
`evaluate_benchmark.py --samplesheet`.

## D38 — The rl100 re-score: `supported` collapses 43 → 5, every remaining gate is identified, and the grading instrument produced a false FAIL for the third time

**Status:** run complete (11h 1m 36s, 150 tasks, exit 0). Evaluator bug fixed; no
threshold changed.

### Tier movement against the Oct 7 baseline (all 11 samples, 82 findings)

| tier | Oct 7 | 2026-10-20 | Δ |
|---|---|---|---|
| `confirmed_in_silico` | 0 | **0** | — |
| `supported` | 43 | **5** | **−38** |
| `attributed_to_neighbour` | 16 | 45 | +29 |
| `screened_only` | 23 | 32 | +9 |

The collapse is the same direction rl150 moved on its Oct 15 re-score (6 → 1) and it
is **specificity improving, not sensitivity lost.** These are single-organism
libraries, so a *Salmonella* call in the *E. coli* library is a cross-attribution
artifact and `attributed_to_neighbour` is the right answer. Tier PPV is now monotone:

| tier | calls | truth present | truth absent | PPV |
|---|---|---|---|---|
| `supported` | 5 | 3 | 2 | **0.600** |
| `attributed_to_neighbour` | 45 | 5 | 40 | 0.111 |
| `screened_only` | 32 | 0 | 32 | 0.000 |

No false confirmations on absent taxa. `confirmed_in_silico` has no calls, so the
question the table exists to answer — does the top tier earn its complexity — still
cannot be answered on real data.

### What holds each of the five `supported` calls

Exactly **one** passed every alignment gate: *S. enterica* in its own library, capped
solely by the missing marker rule block. That is D37's diagnosis confirmed on the
full run rather than on five samples.

The other four are all *L. monocytogenes*, and all four fail the **same** gate:

```
listeria_monocytogenes    rpm 620881  breadth 0.8751  occupancy 0.9525  identity 0.9643
zymo_log                  rpm 537533  breadth 0.8836  occupancy 0.9528  identity 0.9658
bacillus_subtilis         rpm   1817  breadth 0.7391  occupancy 0.9375  identity 0.9615
saccharomyces_cerevisiae  rpm    356  breadth 0.2092  occupancy 0.7552  identity 0.9646
```

`min_read_identity` for illumina/default is **0.98**. A genuine *L. monocytogenes*
isolate at 620,881 RPM with 95% window occupancy misses by 1.6 points.

**This is D35's ONT finding repeating on Illumina, and the floor is again NOT being
lowered.** On ONT every call was capped by the identity floor and the reason was
R9.4.1 chemistry; here the cause is different — Illumina substitution error is
0.1–1%, so 96.4% is **strain divergence from the pack reference**, not sequencing
error. That is a real and expected property of a species reference pack, and it is
the right thing to measure. But 0.98 was chosen as a *species-discrimination* floor
(Wright et al.), and lowering it is exactly the trade this project has refused twice:
it would weaken the gate that separates *E. coli* from *Shigella* and
*F. tularensis* from *F. novicida* in order to admit one organism whose reference is
simply a different lineage. **The remedy is a better reference, or a per-taxon
override, not a lower global floor.** Recorded as the next calibration question.

### The two false positives are probably not false

The two truth-absent `supported` calls are *Listeria* in the *B. subtilis* and
*S. cerevisiae* isolate libraries. They are graded false positives, and the evidence
says they are **real Listeria DNA**: identity 0.9615 and 0.9646 against 0.9643 in
Listeria's own library. Spurious mapping of *B. subtilis* reads onto a *Listeria*
reference would not land within 0.3 points of the real thing, and breadth 0.7391 at
77,052 reads is 74% of the genome, not a trace. Library cross-contamination (index
hopping is documented for this dataset) is the parsimonious reading.

**The benchmark cannot tell that apart from a false positive, because this run has no
negative control** — `samplesheet_illumina_rl100.csv` declares no `role=control` row,
and the report says so. Until it has one, `supported` PPV of 0.600 is a floor, not a
measurement. Adding a control to the real-data samplesheets is now the highest-value
cheap change available.

### The evaluator failed a correct exoneration (third occurrence)

First grading run reported:

```
FAIL zymo_log/Bacillus anthracis max_tier: expected <= screened_only,
                                 observed attributed_to_neighbour
```

*B. anthracis* is absent from the log community by construction, *B. subtilis* IS
present, and they are near neighbours in the pack — so competitive mapping attributed
the reads to *B. subtilis* and said so. **That is the machinery working as designed,
reported as a specificity failure.**

The cause is a structural inconsistency inside `evaluate_benchmark.py`, not a
judgement call about the scenario. `TIER_ORDER` ranks tiers by *how much presence they
assert*, and four tiers that assert none are interleaved with ones that do:

| tier | outcome | rank |
|---|---|---|
| `screened_only` | NEGATIVE | 2 |
| `attributed_to_neighbour` | NEGATIVE | 3 |
| `group_lead` | LEAD | 4 |
| `indistinguishable_complex` | **POSITIVE** | 5 |
| `marker_veto` | NEGATIVE | **6** |

`marker_veto` — the strongest exoneration the pipeline can produce — outranks a
POSITIVE tier. So a rank-only `max_tier` comparison fails the pipeline's best
negatives while passing a weaker positive, and it **contradicts the specificity rate
computed in the same file**, which counts only `PRESENCE_ASSERTING`. The scenario's
own comment already said the narrower thing: *"a confirmed call for any of these is a
false positive."*

`max_tier` is now breached only by a tier that is both above the cap **and**
presence-asserting, and a pass granted on that basis says so rather than going
silently green:

```
PASS  above the cap by rank but asserts no presence: attributed_to_neighbour
      - counted as a correct exoneration, not a false positive
```

Ten new assertions inside `selftest:evaluate_benchmark` pin both halves: a neighbour
attribution, a marker veto and a group lead each PASS a specificity cap and each say
why, while `supported`, `confirmed_in_silico`, `indistinguishable_complex` and
`vaccine_strain` all still FAIL one. No scenario expectation was relaxed.

**Three for three now:** D35 scored a sample that never ran as a miss, D37 printed a
remedy that could not parse, D38 failed a correct exoneration. The instrument that
measures the pipeline has produced more wrong verdicts than the pipeline has.

### `--self-test` was unreachable by hand

`evaluate_benchmark.py --self-test` died on *"the following arguments are required:
--scenarios, --manifest, --findings, --out-dir"*, because argparse validates required
arguments before any code runs — so the suite had been passing four dummy `"x"` values
to get past it. Identical defect and identical fix to `mask_pack_refs.py` in D34,
where `--pack` blocked a self-test that did not exist yet. A self-test nobody can
invoke by hand is one that stops being run.

### The report was regenerated, and nothing moved

`PARSE_SCREEN` completed before D36 was found, so this run's summaries carried the
four false `ARM_EMPTY_REPORT` flags. Regenerated by direct invocation of
`parse_screen.py` and `make_report.py`: **4 flags → 0**, and the "a stage-1 classifier
was enabled and produced nothing" banner is gone, leaving only the true one, *"No
negative control in this run."* `biothreatfinder_summary.tsv` tier columns are
byte-identical before and after — the screen summary reaches `MAKE_REPORT` and never
`SCORE_EVIDENCE` (D33), so it cannot move a call. `score_evidence/findings.tsv` was
not touched. `runs/zymo_rl100_rescore/report/REGENERATED.txt` records this in the run
directory so the regeneration is not invisible to anyone reading the results later.

`-resume` was deliberately not used: changing `bin/` re-hashes the tasks and would
have redone ~11 h of `MARKER_SEARCH_AA` to fix a banner.

## D39 — Declaring a negative control made the report assert the control was clean while reporting that same control as a detection

D38 closed by naming the one cheap change left: *"no real-data samplesheet declares a
negative control, so `supported` PPV of 0.600 is a floor, not a measurement."* Doing
exactly that, on the inputs of a finished run, produced a report that was worse than
the one with no control at all.

### The divergence

`score_evidence.py` asked "which samples are controls" in two places that answered
differently:

| consumer | source of `role` | what it decides |
|---|---|---|
| `control_rpm` | the `role` column on coverage rows | whether anything is **actually subtracted** |
| `n_controls` | coverage rows, candidate rows **and `--samplesheet`** | the QC sentence on every finding |

Inside a pipeline run the two agree, because both descend from the same Nextflow
`meta` — `main.nf:297` parses `role`, `detect.nf` passes `--role` to the scorer, and
the value is written onto every output row. They diverge on the **re-score path**:
published coverage TSVs from a finished run, plus a samplesheet. That is not an exotic
invocation — it is how `zymo_rl100_rescore`, `zymo_rl150_rescore` and the D36 report
regeneration were all produced, and it is the only way to add a control declaration to
a run without paying for it again.

Declaring a control that way set `n_controls` to 1, which silenced *"NO NEGATIVE
CONTROL IN RUN — contamination cannot be excluded"*, and left `control_rpm` empty,
which selected the branch written for the opposite case:

```
1 negative control(s) in run, all clean - no taxon raised in
any control, so no control subtraction applied
```

Measured on `zymo_rl100` with the yeast library declared a control: that sentence went
onto **all 82 findings**, the tier column came back **byte-identical** to the run with
no control at all, and the same file reported the declared control as a `supported`
*Listeria monocytogenes* detection at **355.54 RPM** — one of **14** watchlist taxa it
carries. Two statements in one report, one asserting that no taxon was raised in any
control and the other naming a notifiable organism in it at an actionable tier.

The missing-control case says *"contamination cannot be excluded."* This said
**excluded**, on the strength of evidence it had never read. The branch at fault was
itself written to fix a false alarm — a spotless NTC reading as a missing one — which
is the recurring shape here: the fix for a warning that fires too often becomes a
claim that is false in a case nobody enumerated.

### The samplesheet is the operator's declaration, so it now drives the subtraction

`--samplesheet` is read **first** and relabels the `role` on every candidate and
coverage row whose sample it names, before `control_rpm` is built. It is authority in
both directions, because a stale `control` role on a re-scored row would otherwise
impose a floor the operator never asked for. Every relabelling is reported on stderr
rather than applied silently, and a sample promoted to control is named individually:
28 rows across 1 sample, on the run above.

### A declared control that never ran must not be credited

"No rows anywhere" is the signature of a spotless NTC **and** of a sample the run never
processed — a sheet left over from another batch, a row whose reads failed QC, a name
that never matched what the pipeline emitted. Candidate and coverage rows cannot tell
those apart, which is why the **read ledger** decides it: a control that was sequenced
and analysed has a `seqkit_stats` row even when it nominates nothing. A declared
control missing from the ledger *and* from every row is not counted, and is named on
stderr. Without `--read-stats` there is no ledger, so the question is left unasked
rather than answered wrongly, and the spotless case keeps its credit.

### What the control was then worth

With the subtraction actually running, on the same 11 libraries and the same alignment
evidence:

| | no control | yeast declared control |
|---|---|---|
| `supported` | 5 | 4 |
| `supported` PPV | 0.600 | **0.750** |
| false positives | 2 | **1** |
| false confirmations | 0 | 0 |
| `suppressed_control` | — | 3 calls, **0** truth-present |
| detection rate, every scenario | unchanged | unchanged |

**Nothing true was lost.** All three suppressed rows were truth-absent, and every
per-scenario detection rate is identical — the control removed a false positive and
cost no sensitivity on this data.

The one surviving false positive is *L. monocytogenes* in the *B. subtilis* library,
and it survives on a margin: **1816.73 / 355.54 = 5.11×** against the 5.0× requirement.
The finding discloses it in words — *"control present at 355.54 RPM, sample 1816.73 RPM
(5.1x)"* — and it also carries *"marker rule satisfied: lmono_hly and lmono_prfA"*.
Listeria's own virulence markers, in a *Bacillus subtilis* library. That is not
spurious mapping; it is **real Listeria DNA**, and the contamination hypothesis from
D38 is now quantitatively supported rather than merely plausible. A single control
turns 9% of the specificity error into a disclosed 5.11× margin.

### What was NOT changed

No real samplesheet in the repo was edited. The yeast library is **not** a
no-template control — it is another organism's library from the same sequencing run,
and declaring it one is a benchmark-design decision with a real cost: it imposes a
run-wide floor at 5× whatever hopped into it, which on different data would suppress a
genuine low-level positive. The measurement above is the evidence for making that
call, not the call itself. `benchmarks/scenarios/real_mock_communities.yaml` would need
the control declared alongside it.

### Confirmed on a second real dataset, with the operator's own control

`fsln/case_1/samplesheet_ctrl.csv` declares `case_1-2` a control and has sat beside
`samplesheet.csv` since 2026-10-15 — but `runs/fsln_case_1` was scored with the plain
sheet, so **control subtraction had never run on real data at all**, and that run's
published report carries *"NO NEGATIVE CONTROL IN RUN"* on every finding while a
control declaration sat next to its inputs. Re-scoring with the operator's sheet:

```
case_1-1 | Corynebacterium diphtheriae   screened_only -> suppressed_control
case_1-1 | Toxoplasma gondii             screened_only -> suppressed_control
case_1-2 | Corynebacterium diphtheriae   screened_only -> control_background
case_1-2 | Toxoplasma gondii             screened_only -> control_background
```

The one actionable call — *Candida auris* at 29,705.9 RPM — is **unchanged**. The
control removed two low-level leads and left the real finding alone, which is what a
control is for. `runs/fsln_case_1` and `runs/fsln_case_3` were NOT regenerated: they
are case data and re-publishing their reports is not a documentation change.

Suite at **251 passed, 0 failed, 2 skipped**. The twelve new `ssrole` checks were run
against the pre-fix scorer first: **8 fail, 4 pass**, and the 4 that pass are the
regression guards that must hold in both directions — no-samplesheet runs, the spotless
control, and the no-ledger benefit of the doubt. A test that passes against the bug it
names is the D37 defect, so the check is now part of adding one.

## D40 — The Salmonella marker block: a real `confirmed_in_silico` on real data, and three ways the marker database was holding the wrong gene

D37 established that the binding constraint on the top tier was no longer a threshold
but **curation**: the alignment gates are reachable on real 2×101 Illumina data, and
*S. enterica* in its own library passes all six. What held it at `supported` was D17a's
cap — 32 of 90 watchlist taxa had no `marker_rules.yaml` block. D37 deliberately did
**not** commit a `no_markers_available: true` declaration for *Salmonella*, because the
organism has well-known targets and the declaration would have been false biology.

This is that block, curated rather than declared.

### The block confirms SPECIES IDENTITY, which makes it the opposite of its neighbours

```yaml
"Salmonella enterica":
    confirm_requires: salm_invA and salm_ttrC and salm_ssaV
    support_requires: any_of(salm_invA, salm_ttrC, salm_ssaV)
```

*E. coli* gates confirmation on virulence loci because species presence is commensal
and uninformative. *V. cholerae* gates on `ctxAB` because non-toxigenic *V. cholerae* is
not notifiable cholera. **Neither is true here.** Any *S. enterica* in a food or clinical
specimen is actionable, so gating on a virulence or serovar gene would suppress a
correct and reportable call. Species identity *is* the actionable unit for the detection
decision, even though the serovar drives the public-health response.

| marker | locus | what it is | level |
|---|---|---|---|
| `salm_invA` | SPI-1 | T3SS component; the canonical Salmonella target since Rahn 1992 | **genus** |
| `salm_ttrC` | ttrRSBCA | tetrathionate respiration, centisome 30; Malorny's 2004 real-time PCR target | **genus** |
| `salm_ssaV` | SPI-2 | SPI-2 T3SS apparatus | **species** |

The point of the third row: **invA and ttrC are both genus-level** — *S. bongori*
carries both — so neither alone, nor the pair, says *S. enterica*. `ssaV` supplies the
species call, because *S. bongori* lacks the SPI-2 T3SS (Fookes 2011; SPI-2 was acquired
by a lineage that already had SPI-1). `ttrC` and `ssaV` both sit in the centisome-30
region and so are **not independent loci**; `invA` at SPI-1 is the independent second
location, giving two chromosomal sites across three markers. That distinction is in the
rationale rather than claimed away.

All three are required, so any one untested holds the finding at `supported` under
three-valued logic — the safe direction, deliberately. `support_requires` is `any_of`
because invA-negative *Salmonella* are documented (99.6% inclusivity over 364 strains),
so one absent marker is weak evidence against the organism.

**No `veto_if`, and not invented** — the same discipline as *Y. pestis*. The watchlist
declares *E. coli* (562) and *Citrobacter freundii* (546) as the near neighbours, and all
three markers are absent from both, so the confirm set already discriminates against
everything the watchlist names. The one organism it must exclude is *S. bongori*, and
that is done by **requiring** `ssaV` rather than by vetoing on something *S. bongori*
carries.

### Two things the literature check changed

Written from memory, then checked, and the checks moved the text:

- **The ttr assay is not an ISO method.** ISO 6579 is the *culture* reference method
  that Malorny's PCR was validated against. The catalogue note says exactly that now.
- **`tviB`/`viaB` was considered as a typhoidal marker and REJECTED.** SPI-7/viaB also
  occurs in Vi-positive Paratyphi C and some Dublin, and Vi-negative Typhi are
  documented in endemic regions — so the line would read as "typhoid" to an analyst
  while being neither necessary nor sufficient. Serovar stays out of scope, as the
  watchlist note already asks: *"Serovar is the actionable unit, not species. Hand off
  to a dedicated serovar caller."* A marker that misleads is worse than an absent one.

### The result: `confirmed_in_silico` on real data, on curated evidence

Built the three markers from NCBI (3 built, 0 failed, 0 retrieval errors), aligned the
real `ERR2935848` library against them with the pipeline's own `MARKER_ALIGN_NT`
command, derived states with `marker_evidence.py`, and re-scored:

```
salm_invA  PRESENT  breadth 1.000 >= 0.8, id 0.987 >= 0.95, 12728 reads, 605x
salm_ttrC  PRESENT  breadth 1.000 >= 0.8, id 0.985 >= 0.95,  5101 reads, 475x
salm_ssaV  PRESENT  breadth 1.000 >= 0.8, id 0.989 >= 0.95, 10479 reads, 499x

tier: confirmed_in_silico
  - all confirmation gates passed: breadth 0.9137 >= 0.0100; occupancy 0.9378 >= 0.1500;
    identity 0.9870 >= 0.9800 (illumina/default); depth_cv 0.42 <= 2.50;
    mean_alignment_length 99 >= 61 (illumina 60% of 101 bp read); mapped_fraction 1.000 >= 0.50
  - marker rule satisfied: salm_invA and salm_ttrC and salm_ssaV
  - NO NEGATIVE CONTROL IN RUN - contamination cannot be excluded
  - CAVEAT: KrakenUniq not run, so no distinct-k-mer evidence was available
```

**Exactly one row moved** across all 82 findings: this one, `supported` →
`confirmed_in_silico`. The call still discloses the missing control and the missing
KrakenUniq arm, which is the behaviour that makes the tier worth having.

D37 reached this tier once before, on a **scratch** `no_markers_available` declaration
that was never committed. This is the first time on curated marker evidence.

### Three ways the marker database was holding the wrong gene

The verification kept disagreeing with itself, and each disagreement was a bug.

**1. `salm_ttrC`'s own aliases fetched a different gene.** The row listed the operon
partners `ttrA;ttrB;ttrRSBCA` as aliases. CDS extraction matches any wanted name, so the
nucleotide tier returned **3 × 3063 bp of `[gene=ttrA]`** while the protein tier
correctly fetched **340 aa of ttrC** — one marker, two tiers, two different genes, and
`min_coverage 0.80` measured against the wrong reference. Fixed by removing the aliases:
the rebuild returns 3 × 1023 bp of `[gene=ttrC]`, and 1023 = 340 × 3 + 3, so the tiers
now agree. A descriptive alias would have been unsafe too — a CDS named `tet` would
claim "tetrathionate reductase subunit C" — so the row carries none.

**2. `gene_name_matches` let a one-letter CDS name claim a word alias.** The matcher
accepted `w.startswith(g)` with no length condition. Measured against `HEAD`:

```
True   [gene=P] claims ['L', 'polymerase']
True   [gene=E] claims ['E2', 'TC-83']
```

`"polymerase".startswith("p")`. So **`nipah_L`'s nucleotide reference is 2130 bp of the
Nipah P (phosphoprotein) CDS**, `CAD92360.1`, while its protein tier holds 8 × 2244 aa
of the real L polymerase. `nipah_L` is in `confirm_requires` for Nipah henipavirus, a
BSL-4 agent. Six markers share the L/polymerase shape — `cchf_L`, `marburg_L`,
`lassa_L`, `junin_L`, `rvfv_L`, `nipah_L` — and `nipah_L` is the one that took the CDS
route, so it is the one that was hit. The direction that was *intended* to be permissive
— a CDS name that EXTENDS a wanted name, `capBCADE` for `capB`, `ORF1ab` for `ORF1a` —
is untouched; the reverse now requires three characters, so `pag` for `pagA` still
works and a single letter cannot.

**3. Nothing preferred the gene actually asked for.** `ypest_pmt1_caf1` declares `caf1`,
the F1 capsular antigen, and its nucleotide reference holds **caf1A, caf1M and caf1R —
the operon's accessory genes — and not caf1**, because each is a legitimate prefix match
and nothing ranked them. It is in *Y. pestis*'s `confirm_requires`. `gene_match_rank`
now returns a rank and CDS extraction keeps only the best one present, so an exact match
wins and the permissiveness survives as a fallback.

### A fourth, in the consumer: the states file's columns depended on sort order

`marker_evidence.py` took its header from `results[0]`. A record only carries
`nt_depth`, `nt_expected_breadth`, `nt_n_seqs` and `nt_best_ref` if it reached the
nucleotide measurement block, so a catalogue whose first row is UNKNOWN produced a
header without them — and `extrasaction="ignore"` then dropped those four from **every**
row, including the ones that had measured them, with no error. Observed directly here: a
three-marker database put `banth_pxo1_pagA` (UNKNOWN) first and the states file lost all
four columns while the three PRESENT rows had depths of 474–605×.

`nt_expected_breadth` is the input to the depth gate — the rule that marker absence is
only called when depth could have shown it. The column that exists to stop a false
ABSENT could vanish silently. Fieldnames now come from the union of every record. Some
keys were already pre-initialised against exactly this, with a comment saying so; that
fix was per-column and could not cover a block adding four at once.

### The invariant that found two of these, now a test

Declared gene versus the `[gene=]` tag on the retrieved CDS, across the whole shipped
marker database, in two strengths:

- **no built reference is a gene the matcher would reject** — judged by the builder's
  own `gene_match_rank`, so the invariant asks the same question the retriever asked.
  `ORF1ab` for a wanted `ORF1a` passes; the Nipah `P` does not.
- **every built reference contains the gene actually declared** — exact match on the
  gene *or* a declared alias. `wbeT` for `wbeO1` and `MP` for `M` are right answers and
  pass; what this catches is a marker holding only its operon's other genes.

115 of the 166 catalogue rows have a CDS gene tag to check (31 took the whole-record
route and have none; 17 have no nucleotide file). Two entries are pending, each with its
cause, and the list is **self-retiring**: the test also asserts that every pending defect
is *still real*, so a rebuild makes it fail and forces the entry out. A pending item
cannot outlive the problem it describes.

### What was NOT done

- **The production marker database was not rebuilt.** `nipah_L` and `ypest_pmt1_caf1`
  are still wrong on disk; the code that produced them is fixed, and realising that
  needs a full rebuild of all 166 markers. `--only-marker` cannot do it safely —
  `concat_files()` rewrites `markers.nt.fna` from the manifest rows while the index
  builders skip when their output exists, which is the hazard recorded on
  `md5_or_none()` and the reason `bmand_18S` once came back ABSENT rather than UNKNOWN.
  Everything above was built into a scratch `--db-root`.
- **No samplesheet or benchmark file was changed**, so D39's control decision is still
  open and this run still carries "NO NEGATIVE CONTROL IN RUN".
- **`d17a` was repointed, not relaxed.** It asserted the no-block cap using
  *Salmonella* as "a taxon which really has no block", and curating *Salmonella* made
  that premise false — so the test went red, correctly. It now uses *Orientia
  tsutsugamushi*, as of 2026-10-20 the only notifiable bacterium left without a block.
  31 of 90 taxa still have none.

## D41 — The marker database rebuild: a broken local library degraded a third of the catalogue, and the first attempt was rolled back

D40 fixed the retrieval code but left the shipped marker database untouched, because
realising the fix needs a full rebuild of all 166 rows and `--only-marker` cannot do it
safely. This is that rebuild. It took two attempts, and the first one had to be undone.

### `--force` was not optional

Without it the resume path SKIPS any row whose status is `built` or `partial` when the
files are present, the md5 matches and no catalogue column has drifted — which is
exactly the state of the two markers the rebuild existed to fix. A plain re-run would
have reported success and changed nothing. `--force` also bypasses the 13 GB Entrez
response cache, so all 166 rows were re-fetched from NCBI.

### The first attempt degraded 33 markers, and the cause was local

It completed, exited non-zero, and the comparison against the pre-rebuild manifest read:

| | before | after attempt 1 |
|---|---|---|
| `built` | 143 | **113** |
| `partial` | 4 | **37** |
| improvements | — | **0** |

Every one of the 33 regressions had the same cause, and it had nothing to do with NCBI
or with the retrieved sequences: `hmmbuild` exited 127 with
`libgsl.so.25: cannot open shared object file`. The profile step failed for every marker
declaring an HMM tier, and a marker missing a declared tier is recorded `partial`. The
sequences were fine; the profiles simply never got built.

**The database was restored from the backup taken before the run** — 143 built, 33
profiles, 163 rows, verified against the saved manifest. The regressed copy was kept
rather than deleted, so the failure can be inspected instead of reconstructed.

Two process lessons, both mine. The PATH was assembled by searching for the first
executable of each name, which found an `hmmbuild` that resolves but cannot run — so
the pre-flight check reported it as present. **Resolving a tool is not the same as
running it**, and the dry run printed the path without ever invoking it. And the second
attempt was launched after validating the replacement on a single profile marker in a
scratch `--db-root` first, which is what the first attempt should have done.

### The second attempt is clean

| | baseline | rebuilt |
|---|---|---|
| `built` | 143 | **146** |
| `partial` | 4 | 4 |
| `failed` | 3 | 3 |
| `motif_todo` | 9 | 9 |
| `curation_required` | 4 | 4 |
| profiles (`has_hmm`) | 33 | 33 |
| status regressions | — | **0** |

The +3 are D40's Salmonella markers, which until now existed only in a scratch
database. Zero profile failures. Every other category identical to baseline.

### Verified against the production database rather than a scratch copy

The real library was realigned against the rebuilt `markers.nt.fna`, states were derived
from the production manifest, and the run was re-scored: all three Salmonella markers
PRESENT at breadth 1.000 and identity 0.985–0.989, and the finding reaches
`confirmed_in_silico` on `marker rule satisfied: salm_invA and salm_ttrC and salm_ssaV`.
D40's verification used a three-marker scratch database; this one uses the database the
pipeline will actually load.

### Both pending gene-tag defects are resolved, and one of them not as predicted

`ypest_pmt1_caf1` now holds the gene it declares, 513 bp, where it had held only the
operon's accessory genes.

`nipah_L` was predicted to stay broken and did not. It no longer holds the mismatched
coding sequence: with the corrected matcher nothing matched at the CDS level, so it fell
back to the whole-record route — and that route is intended, since
`NT_MAX_DIRECT_LEN = 50_000`. Compared against its five sibling markers, it now has the
same strategy, the same coverage floor, the same evidence tiers and the same retrieval
route. It is consistent with the pack's design instead of being the one outlier holding
a different gene, so the defect D40 identified is gone.

The self-retiring pending list did its job: it went red the moment the rebuild made both
entries stale, and each was checked against the rebuilt database before being removed
rather than deleted because the test complained. `GENE_TAG_MUST_FIX` is now empty.

### Still open, and NOT introduced here

16 markers that declare a gene hold a whole viral genome or segment as their nucleotide
reference, because the whole-record route accepts anything under 50 kb. Their coverage
floors are 0.50, so for those rows the nucleotide tier asks "is half of this genome
covered" rather than "is this gene present", with the protein and profile tiers carrying
the specificity. That is pre-existing, unchanged by this rebuild, and a question about
marker semantics rather than a bug to fix quietly. Recorded rather than touched.

Suite **286 passed, 0 failed**.

## D42 — A second reference for one species: identity 0.9641 → 0.9947, 48% more reads recovered, and the floor stays at 0.98

D38 measured a genuine *L. monocytogenes* isolate at identity 0.9643 against the 0.98
Illumina floor and refused to lower it, because 0.98 is a species-discrimination floor
protecting pairs the pipeline must keep apart. The remedy named there was "a better
reference or a per-taxon override". The question asked here was whether a species can
carry more than one reference. It can, the pack already does it for viruses, and it
turns out to be the right answer.

### The pack already encodes this reasoning - for viruses only

`build_pack_refs.py` takes `--max-per-taxon` (default **1**) and
`--max-per-taxon-virus` (default **3**), and the help text gives the rationale:
*"Viruses need more: intra-species diversity is high enough that one RefSeq genome can
lose reads from a divergent strain."* 74 of the pack's taxids carry two or three
references, every one of them a virus; two notifiable viruses carry three target
assemblies each. **No bacterium carried more than one.** The Listeria measurement is
evidence that the same reasoning applies to some bacteria, so this is the existing
design extended, not a new mechanism.

### The measurement, same reads both ways

The isolate library aligned against the one-reference bundle and against the
two-reference bundle, same tool, same MAPQ, `coverage_stats.py` for both:

| metric | 1 reference | 2 references |
|---|---|---|
| `mean_identity` | 0.9641 | **0.9947** |
| `best_reference_breadth` | 0.8768 | **0.9822** |
| `best_reference_occupancy` | 0.9525 | **0.9914** |
| reads mapped | 12,826,145 | **19,011,128** |
| raw `breadth` | 0.8768 | 0.5444 (diluted) |
| raw `window_occupancy` | 0.9525 | 0.7436 (diluted) |

The one-reference figure of 0.9641 reproduces the production run's 0.9643, which is the
control that makes the comparison valid: the production run aligned only the ~1.46 M
screen-assigned reads while this used the whole library, and identity is unchanged by
that difference.

**Three things improved, not one.** Identity clears the floor by a wide margin. The best
reference is *better* covered, 0.877 → 0.982, because the correct-lineage genome is
covered end to end rather than patchily. And **6.2 million more reads map at all, a 48%
gain** - they previously could not align well enough to the wrong lineage. That is
recovered sensitivity, not a threshold workaround. The reads overwhelmingly prefer the
newly added assembly, so the isolate genuinely is the lineage the pack lacked.

**The 0.98 floor is untouched.** No per-taxon override was added, and none is needed.

### Dilution is a property of having several references, not of rank

`breadth` and `window_occupancy` are covered-over-SUMMED across every reference filed
under a taxid, so the second genome halves both - visible above as 0.877 → 0.544.
`coverage_stats.py` already emitted `best_reference_breadth` and
`best_reference_occupancy` on **every** row for exactly this, but `score_evidence.py`
substituted them for supra-specific entries only. Correct as far as it went, and wrong
about why: the dilution follows from having more than one reference, whatever the rank.

The substitution now applies to any multi-reference taxon, with one condition.
`best_reference_*` is measured per BAM reference, i.e. **per contig** - so for a draft
assembly the "best reference" is the best-covered contig, and substituting it would let
a well-covered 40 kb contig stand in for a 5 Mb genome. It is exact only when every
assembly under the taxid is a single sequence, which `load_ref_bundles()` answers from
the manifest's `n_contigs` without re-running coverage. When that fails the diluted
numbers are kept and the row says it may be UNDER-called - under-calling is the safe
direction, but a reader must not have to work out why the floors were applied to a
denominator including genomes the sample never had.

### Two bugs, both the shape this project keeps finding

**`--max-per-taxon N` could not produce N distinct assemblies.** The candidate dedupe
compared **filenames**: a reused genome sits under the plain `<taxid>_<name>.fna.gz`
while an extra is written as `<taxid>_<name>.<accession>.fna.gz`, so re-selecting the
accession already on disk produced a different filename, passed the check, and was
downloaded again. Measured: four taxa ended with two manifest rows naming **one**
accession and two byte-identical copies of one genome, while the summary reported
*"targets fully built: 1 / 1"* and *"genomes on disk: 8"*. The flag was unusable for any
taxon that already had one reference - which is every taxon in a built pack. Now deduped
on accession, and the rerun fetched a genuinely distinct assembly.

**Adding a reference did not invalidate the competitive bundle.** The FASTA and its
minimap2 indexes are skipped whenever the output merely *exists*, so the manifest said
two references while the aligner still used one - a disagreement invisible from either
side. Existence is not freshness. Rebuilt explicitly here; the general fix is a
staleness check against the inputs, and it is not done.

### What is now pinned

Three invariants over the shipped pack, each verified to FAIL against a reproduction of
the broken state: a `(taxid, accession)` pair appears once; a taxon's row count equals
its distinct accession count; and no competitive bundle is older than the genomes it is
built from. Plus eight checks on the geometry substitution, six of which fail against
the pre-change scorer.

The masking invariant from D34 did its job unprompted: the suite went red with *"4
unmasked"* the moment the new references landed, and stayed red until they were masked.
375 references, all with a masking result.

Suite **297 passed, 0 failed**.

### Not done

The three near-neighbour taxa in that bundle also gained a second assembly, because
`--max-per-taxon` applies to neighbours too. That improves competitive absorption and
costs nothing measured here, but it was a side effect rather than a decision. No other
bacterium was touched: the question of which bacteria deserve a second reference is now
answerable by measurement, and has been answered for one.

---

## D43 — The worst-case reference measurement: how far the shipped references reach across each target species, measured without a specimen

D42 gave one species a second reference and showed identity 0.9641 → 0.9947 on real
reads. The measurement was sound and the way it was reached was not: **the species was
chosen because a sample complained, and the improvement was scored on the same reads
that motivated it.** No isolate was held out, and the other target bacteria were never
examined, because no specimen had yet accused them. For a surveillance pack that is
backwards. A true sample is an isolate nobody has seen.

So this asks a question that needs no specimen at all: *for every assembly of this
species that is publicly known, would reads from that isolate clear the gates against
the reference set we ship?* The worst answer over all known assemblies is the pack's
guarantee. It is a property of the pack, not of a run.

### Scope and method

Frames were taken as the RefSeq complete/chromosome subtree of each of the 26 bacterial
target taxids — **13,582 assemblies**, all of them downloaded, so the result is a true
minimum rather than an estimate from a sample. Each was sketched (mash, k=21, s=10,000)
and measured against the **masked** shipped FASTA, because that is what reads are
actually aligned to.

`scripts/measure_reference_span.py` carries the whole chain — `frame`, `fetch`, `span`,
`cover`, `verdict` — so it is repeatable rather than a one-off. Its `--self-test` is in
the suite.

### Why the verdict is not an ANI number

ANI is cheap enough to run over everything, so it is used **only to choose which genomes
to test** and to count how many sit near the edge. It cannot carry the verdict, because
the ANI at which the pipeline's `mean_identity` crosses the 0.98 floor is
**species-dependent**, and not by a little:

| species | `mean_identity` crosses 0.98 at mash-ANI |
|---|---|
| *L. monocytogenes* | ~0.974 |
| *E. coli* | ~0.973 |
| *L. pneumophila* | ~0.981 |

Read identity is measured over the regions that still align, so how far it floats above
genome-wide ANI depends on whether a species' divergence is uniform drift or a few
replaced loci. A single ANI cut-off across the pack would therefore under-call some
species and over-call others. The `verdict` stage instead simulates reads from the
assemblies ANI nominated as worst and pushes them through **the pipeline's own** aligner,
competitive bundle and `coverage_stats.py`. The deciding number is the one the pipeline
would itself produce. Calibration panels were measured, not assumed.

Two instrument controls passed before any of it was trusted: each shipped reference
scores ANI 1.0000 against itself, and a genome simulated against itself returns
`mean_identity` 0.9980 — not 1.0, because of the read-error model. The ceiling is
0.998 and the floor is 0.98, so there are under two points of headroom.

### What the measurement found

**Most of the pack is fine, and that is a real result.** 17 of the 26 species need
no addition: their most divergent known assembly already sits within reach of the single
shipped reference. Several are clonal enough that one reference is genuinely sufficient.
The D42 treatment is **not** broadly necessary, which is worth knowing before applying it
anywhere else.

**Nine species do not span their own diversity**, and the shape differs:

- For most, divergent members lose identity and breadth gradually — the finding settles
  below `confirmed_in_silico` and reads like a cautious call on a weak signal rather than
  like a reference gap.
- For one, the degradation is not gradual. Nine assemblies share **under 100 of 10,000
  hashes** with the shipped reference, which is not a low ANI but **no measurable
  similarity**. Simulated reads from those isolates return breadth **0.002** — below even
  the permissive 0.01 floor — from a library of ~360,000 reads. That is not a downgrade.
  It is indistinguishable from a clean negative.

For the D42 species itself the measurement also settles the open question: with the
single original reference, **47.1% of 853 known assemblies** would have failed the
identity gate. The sample D38 ran into was not an unlucky outlier; it was a coin flip.
With both references that falls to **3.6%**, and those remaining sit in the two
genuinely divergent lineages. Two more references would close it. The second reference
therefore did real work — and the measurement, not the sample, is what shows it.

### Where this intersects the tier ladder

`score_evidence.py` requires the coverage-geometry gates to **pass** for a finding to be
actionable (the reasoning is recorded at the `actionable_requires` branch), and marker
rules only ever hold a finding down — they never lift one. So where a shipped reference
cannot recruit reads from a divergent member of its own species, the row is
**not actionable regardless of marker evidence**. For the species above whose watchlist
note already says the actionable signal is a gene cluster rather than the chromosome,
that design intent is sound, but it means the marker route is a **single point of
failure** for those isolates: there is no corroborating reference evidence to fall back
on, so a marker gap and a true negative look the same.

### The selection defect this also exposes

`build_pack_refs.py` ranks candidates on NCBI metadata alone — `refseq_category`, then
assembly level, then GCF over GCA, then size (`Candidate.sort_key`). **Nothing in that
key knows what a lineage is.** `--max-per-taxon 2` can therefore return two assemblies
from the same clade and leave the other half of a species uncovered. D42's second
reference happened to land in the other lineage; that was luck in the pick, not a
property of the picker. The `cover` stage replaces the luck with a stated criterion:
greedy set cover over the species' own diversity, **seeded with what the pack already
ships**, so its output is a list of additions and never churns a reference that earlier
measurements were made against. Ties go to single-sequence complete assemblies, which
D42's geometry substitution needs.

Across the nine species the cover asks for **32 further references**. One species
accounts for nine of them because it is polyphyletic; its watchlist note already says so.

### What is pinned

`--self-test` covers the four things that could produce a confident-looking table rather
than an error: the plain shipped filename may be claimed only by a taxid's **first**
reference (handing it to a later accession would measure the wrong genome and report it
under the wrong name — silently, and in the flattering direction); the cover is seeded,
reports additions, and prefers single-sequence assemblies on ties; a near-zero
shared-hash count is counted as **unmeasurable** rather than reported as a precise low
ANI; and `require_tool()` **runs** every binary before use. That last one is D41's lesson
encoded — that rebuild degraded a third of the marker catalogue because a tool resolved
to a path and then exited 127 on a missing library, and the dry run printed the path
without ever invoking it. The resolver guard was verified to fail against its own first
version.

### Not done

No reference was added. This decision is the measurement and the recommendation; acting
on the 32 is a separate change, species by species, and three of the nine species are
ones whose watchlist notes say species-level alignment is not the intended route anyway —
for those, the right response may be to strengthen the marker route rather than the
reference set.

The report is `docs/reference_span.tsv`. It cannot be regenerated in CI: the frames are
59 GB of downloaded assemblies, so the committed TSV is the evidence and the script is
how to redo it. A CI check against the committed worst-case figures is not built.

One further observation, flagged rather than chased: for a species group whose members
are near-identical, simulated reads from one member spread across the sibling references
in the same bundle, so breadth lands near 0.5 at an ANI of 0.999. That is competitive
mapping working as designed, but it means breadth can be diluted by siblings rather than
by divergence, which is a different failure mode from the one measured here.
