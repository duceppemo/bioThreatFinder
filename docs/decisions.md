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

**Falsifiable prediction, now sharper than the one D12 first recorded.** Reads scale
linearly with spike level (2,114 true reads at 100 rpm), so 3,840 is reached at **~180 rpm**.
Detection should switch on between 100 and 1000 rpm and specifically near 180 - which is
exactly why the measured LoD came out as 1000: the series steps 100 -> 1000 and brackets
the real threshold. A finer grid should find the cliff at ~180 rpm, NOT at 1000. If it is
found anywhere else, this arithmetic is wrong.

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

## D13 — Masking is detected at 200 bp / 95% identity but reads map at 101 bp; the two scales do not meet

**Finding:** `[MEASURED]` 2026-10-03, from the negative-control bundles of
`runs/lod_v2_kraken2` (no watchlist target spiked; reads carry their source accession).
Reads at MAPQ>=10 on each TARGET reference, by true origin:

| taxid | taxon | reads on target | dominant true origin |
|---|---|---|---|
| 5811 | *Toxoplasma gondii* | 24,101 (7,645 in top 3) | **human** NC_060925.1 / 060933.1 / 060943.1 |
| 487 | *Neisseria meningitidis* | 317 | *H. haemolyticus* (275), *E. coli* (42) |
| 1392 | *Bacillus anthracis* | 66 | *L. innocua* (64) |
| 1639 | *Listeria monocytogenes* | 45 | *B. thuringiensis* (45) |
| 622 | *Shigella dysenteriae* | 61 | *E. coli* (60) |
| 28901 | *Salmonella enterica* | 23 | *H. haemolyticus* (20) |

The *T. gondii* reference is a **human-read sink**: thousands of host reads that survived
depletion map to it at MAPQ>=10, three orders of magnitude more than any bacterial
reference recruits. Its neighbour side holds 33 reads, so `neighbour_ratio` is **730** and
the neighbour gate offers no protection whatsoever. Breadth, occupancy and identity are
what actually hold the line - the specificity table shows 0.000 false positives - so this
is a redundancy that is doing nothing, not an active false positive.

**Why masking does not catch it.** `mask_pack_refs.py` runs at `min_hit_len=200` and
`min_identity=0.95`. A 101 bp read needs neither: it maps at MAPQ>=10 on far shorter and
far less identical stretches than a 200 bp 95%-identity alignment. The result pack-wide:

* *T. gondii*: **1,982 bp** host-masked on a 65.7 Mb genome
* **2 of 359 genomes** have ANY host masking at all, 2,737 bp in total
* `crossgroup` is **0 bp** for every genome checked, because it compares only across
  distant groups by design

So the masking step cannot see read-scale recruitment. This is a scale mismatch, not a
missing source.

**Consequence for D12 option 1.** Masking within-family recruitment was listed there as the
lever that would fix numerator and denominator together. At current parameters it would
not: the *E. coli* pedestal on *Y. enterocolitica* sits on a reference with 571 bp masked
and 0 bp crossgroup. Making masking catch it means read-scale parameters (~100 bp, lower
identity), which on closely related genomes would mask large fractions of the reference and
cost sensitivity in exactly the silent way masking exists to prevent. That trade-off has to
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

**The reads are 100% human.** Read attribution on the bundle BAM: of 8,380 reads on the
*T. gondii* reference at MAPQ>=10, **8,380 come from T2T human chromosomes** (NC_0609*).
Not a majority - all of them. 44% of a 65 Mb parasite genome is "covered" at 96.5% window
occupancy and 98.1% identity by host reads that survived depletion.

**It is systematic, not a *T. gondii* quirk.** Measured recruitment per million human reads
(`scripts/measure_recruitment.py`):

| reference | records per 10^6 human reads | % of human reads |
|---|---|---|
| *Toxoplasma gondii* | 38,677 | **3.87%** |
| *Plasmodium knowlesi* | 22,170 | 2.22% |
| *Entamoeba dispar* | 4,080 | 0.41% |
| *Plasmodium vivax* | 2,333 | 0.23% |
| *Trichinella pseudospiralis* | 1,650 | 0.17% |
| *Zika virus* | 1,082 | 0.11% |

⚠️ **These figures were first published at 2x these values** (*T. gondii* 77,354 / 7.74%).
`measure_recruitment.py` divided SAM records, which count each mate separately, by read
PAIRS - so it reported per-pair and called it per-read. Corrected 2026-10-03 after a review
pass; the tool now emits `rpm_source_reads` against `source_reads = 2 x pairs`, matching
`bin/coverage_stats.py:305`, whose denominator sums FASTQ records across both mates. The
RANKING and the mechanism are unaffected, because both were direct record counts.

Eukaryotic references absorb host reads at percent-level rates. The *Zika* row deserves
separate attention: 0.22% of human reads landing on a 10.8 kb RNA virus is not plausible
homology and suggests human or vector sequence INSIDE that reference - the Wright et al.
mechanism, which would make it a contaminated reference rather than a homology artefact.
Flagged, not concluded.

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

**What is actually holding the line, and it is not what the design says.** Back-calculating
from the corrected coefficient: 8,380 recruited records / 0.0387 => **~217,000 surviving
human reads**. Recruitment scales linearly with whatever survives, so **host-depletion
efficiency is the real control on eukaryotic false positives**, and until 2026-10-03 nothing
in the pipeline recorded it (`ch_host_stats` was a dead channel; now written to
`<outdir>/qc/host_depletion.tsv`). It still GATES nothing.

⚠️ A matching prediction was claimed and does not hold. Commit `276a3f1` reported the
measured coefficient predicting the observed *Y. pestis*-bundle pedestal "within 8%". That
compared a per-PAIR coefficient against a per-RECORD count. With consistent units the
prediction is 1,266 records against 2,664 observed - **2.1x out**, not 8%. The direct
observations in this entry stand; the predictive arithmetic does NOT, and the most likely
residual cause is that `unique_reads` for the background organism underestimates the
organism's true record count. Resolving it needs the read-attribution denominator, not the
MAPQ-filtered one.

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
