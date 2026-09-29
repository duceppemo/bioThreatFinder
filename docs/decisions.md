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

**Nine** distinct bugs in this project have shared one shape: *a mismatch that produces an
empty or degraded result indistinguishable from a clean negative.* Rotted taxid · GTDB
namespace · refs base path · crossgroup subset · duplicated bundle sequence · stale
taxonomy · kraken2 2.0.8 lacking the minimizer flag · marker manifest truncation · refs
manifest truncation.

So the first question for any new component is **"could this silently produce nothing?"**
— and the answer must be a guard that fails loudly, plus a check against the written
artifact rather than an assertion in code. Two of the nine were found only by re-reading
what was written; one fix was itself vacuous and passed its own logic.
