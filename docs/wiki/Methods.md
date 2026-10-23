# Methods

## The workflow

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/duceppemo/bioThreatFinder/main/assets/workflow/bioThreatFinder-workflow-dark.svg">
  <img alt="bioThreatFinder workflow" src="https://raw.githubusercontent.com/duceppemo/bioThreatFinder/main/assets/workflow/bioThreatFinder-workflow-light.svg" width="100%">
</picture>

Boxes are named as they appear in the Nextflow trace. `[opt]` is off by default, and
`x candidate` marks the fan-out: the competitive alignment runs once per (sample ×
watchlist candidate), while the marker tier runs once per sample off the host-depleted
reads. Regenerate the diagram with `python3 scripts/make_workflow_svg.py`.

<details>
<summary>The same diagram as text</summary>

```text
  samplesheet.csv ──┬── illumina (paired) ──► FASTP ──┐
                    └── ont (single)      ──► CHOPPER ─┤
                                                       ▼
                                      MINIMAP2_HOST_DEPLETE          pack.yaml
                                      host reads removed                 │
                                   ┌───────────┴───────────┐             ▼
                                   │                       │     RESOLVE_WATCHLIST
                      qc/host_depletion.tsv          SEQKIT_STATS   name ──► taxid
                      how much host got through           │         + descendants
                                                          │               │
  ═════════ STAGE 1 ─ what MIGHT be here (permissive, any-arm OR) ════════│════════
                                                          │               │
                    ┌──────────────┬──────────────┬───────┴──────┐        │
                    ▼              ▼              ▼              ▼        │
                 KRAKEN2        SYLPH_PROFILE  METABULI      KRAKENUNIQ   │
                 + BRACKEN      ANI-based      DNA + AA      distinct     │
                 (illumina)                    [opt]         k-mers [opt] │
                    └──────────────┴──────────────┴──────────────┘        │
                                          │                              │
                                          ▼                              │
                                    PARSE_SCREEN ◄───────────────────────┘
                                    candidates.tsv per sample
                                          │
  ═════════ STAGE 2 ─ does the ALIGNMENT hold up? (strict, this decides) ══════════
                                          │
                 ▼
    MINIMAP2_COMPETITIVE [per candidate]            2b. MARKER TIER
                                                    ▲ per sample,
                                                    │ not per candidate
       clean reads ─────────────────────────────────┘
    target + near neighbours in ONE index           MARKER_ALIGN_NT  (nucleotide)
                 │                                  MARKER_SEARCH_AA (protein)
                 │                                          │ hit reads
                 │                                          ▼
                 │                                  MARKER_SEARCH_HMM (profile)
                 ▼                                          │
          COVERAGE_STATS                                    ▼
          breadth · window occupancy               MARKER_EVIDENCE
          depth_cv · identity                      PRESENT / ABSENT / UNKNOWN
          neighbour_ratio                          per virulence locus
          masked regions excluded
                 │                                            │
                 └────────────────────┬───────────────────────┘
                                      ▼
                               SCORE_EVIDENCE
                               gates + marker rules + control subtraction
                                      │
                                      ▼
                                 findings.tsv
                                      │
                                      ▼
                                 MAKE_REPORT
```

</details>

## Why two stages

Stage-1 classifiers answer "what is the composition of this sample". The question here is
"is organism X present", which has a far harsher false-positive cost. A handful of k-mers
shared with a common environmental relative is enough to place a high-consequence taxon in
a Kraken2 report.

Stage 2 maps each candidate's reads against a bundle containing the **target plus its near
neighbours**, then asks four things k-mer counting cannot:

- **Breadth** — what fraction of the reference is covered at all? A real organism at low
  abundance spreads thin across the whole genome; an artifact stacks on one conserved gene.
- **Window occupancy** — breadth's robust sibling. At very low depth breadth is small
  either way, but 40 reads across 40 windows is an organism and 40 reads in 2 windows is a
  repeat.
- **Competitive attribution** — do these reads actually prefer the target, or its harmless
  relative? `neighbour_ratio` below ~1 means the relative.
- **Control subtraction** — is this above what the negative control shows?

## Evidence tiers

Every finding carries a tier **and** a `reasons` string naming the gate that decided it.
Nothing is a bare verdict.

| Tier | Meaning |
|---|---|
| `confirmed_in_silico` | Met every threshold. **Still requires wet-lab confirmation.** |
| `supported` | Real alignment evidence, below the confirmation bar |
| `no_reference` | On the watchlist and raised by the screen, but no reference genome available — **a database gap, not a negative result** |
| `suppressed_control` | Comparable level in a negative control; contamination until proven otherwise |
| `indistinguishable_complex` | One of a named pair is present and alignment cannot say which — *Shigella* within *E. coli*, *C. deneoformans* within *C. neoformans*. Both references cover well, so a one-sided attribution would be a confident wrong answer. Reachable only with `--contest_attribution`, which **defaults off**: the uniformity filter that separates a cross-mapping pedestal from a genuine sibling turned out to be abundance-driven and is parked pending recalibration |
| `group_lead` | A **group-rank** watchlist entry (family or genus) cleared its gates. A member of the group is present; which member is not a question this evidence can answer, so it is not a species call |
| `attributed_to_neighbour` | Reads align better to a relative |
| `control_background` | What the control itself contains (QC, not a finding) |
| `screened_only` | A stage-1 lead that failed stage 2 — the expected fate of classifier noise |

## Stage 1: the classifier arms

One classifier is not enough for an "is organism X present" question, and the arms fail in
different directions. **Two run by default; two are opt-in.**

| Arm | Default | What it contributes |
|---|---|---|
| **Kraken2 + Bracken** *(primary)* | on | Distinct **minimizer** count and duplicity via `--report-minimizer-data`, plus Bracken abundance |
| sylph | on | ANI/containment; very low false-positive rate for bacteria and fungi |
| Metabuli | **off** | Joint DNA + amino-acid k-mers, for divergent agents. Measured, it has changed **no** finding on any scenario run |
| KrakenUniq | **off** | Distinct **k-mer** cardinality. A second opinion, not a dependency |

**Combining rule: any-classifier OR.** Any arm raises a candidate; stage 2 discriminates.
Requiring consensus would discard the sensitivity that running several arms buys. What
consensus *does* do is get recorded — `n_classifiers` and the arm list travel into scoring
and into the report, because a taxon seen by three independent methods is much stronger
evidence than one seen by a single k-mer counter. K-mer quality is likewise **flagged,
never enforced** at stage 1: a flag that reaches scoring is more useful than a candidate
that silently vanished.

### Why Kraken2 is primary and KrakenUniq is not

An earlier version of this design made KrakenUniq primary, on the grounds that "Kraken2
only counts reads". That was wrong, and reading the source settled it: Kraken2's
`--report-minimizer-data` columns are *literally KrakenUniq's HyperLogLog implementation* —
same author, same `hyperloglogplus.h`, same default precision p=12 — counting distinct
**minimizers** instead of distinct **k-mers**. The only published comparison (Oskolkov
2026) finds **comparable F1** once the threshold is shifted by the 2–3× that minimizer
subsampling implies: ~1,000 unique k-mers for KrakenUniq against ~2,000–3,000 distinct
minimizers for Kraken2. KrakenUniq has had no commit since December 2023.

So the distinct-k-mer signal is kept and KrakenUniq's index is not. If you want it as a
second opinion, **download a prebuilt index rather than building one** — a build costs
~1400 GiB transient and needs Jellyfish 1.x. A Kraken2 index **cannot** be reused:
KrakenUniq descends from Kraken 1 and uses its format, while Kraken2's hashed-minimizer
structure is lossy with respect to the k-mers KrakenUniq must count.

Three caveats on Kraken2's columns, verified from source rather than inferred: columns 4–5
are **clade-aggregated**, so threshold at species rows and not genus rows; minimizers are
keyed by the minimizer's database LCA, not the read's final call; and upstream issue #445
(column 5 exceeding `kraken2-inspect`'s count) is still open, so treat column 5 as a
monotone evidence score to calibrate empirically, not as a literal breadth measurement.

The current flag threshold is `min_distinct_kmers = 2000`, taken from Oskolkov's measured
Kraken2 optimum rather than reasoned; it still needs recalibrating on your own panel and
depth.

### Two guards worth knowing about

**Taxonomy namespace.** Metabuli publishes GTDB-based indexes as well as NCBI-based ones,
and **GTDB taxids are not NCBI taxids**. Point the pipeline at a GTDB index and every
watchlist lookup misses — producing not an error but an empty result, indistinguishable
from a clean negative across the entire screen. So the pipeline reads
`BTF_BUILD_MANIFEST.txt` and **refuses to run** on a non-NCBI index rather than warning.
`parse_screen.py` detects the same mismatch from the data side, in case an index arrives
without a manifest.

**HyperLogLog cardinalities are not additive.** Summing the `kmers` column across
descendant nodes double-counts shared k-mers and inflates the number, sometimes badly. So
distinct-k-mer evidence uses the watchlist node's own clade estimate, falling back to the
*maximum* over descendants — never the sum. Read counts, being plain integers over disjoint
assignments, are summed normally.

## The marker tier (stage 2b)

For a large share of the watchlist, **the chromosome is not the actionable signal**, so
alignment metrics alone cannot produce a defensible call:

| Organism | What actually matters | Why breadth is not enough |
|---|---|---|
| *B. anthracis* | pXO1 toxin + pXO2 capsule | Chromosome is near-identical across the B. cereus group |
| *Y. pestis* | pPCP1 / pMT1 plasmids | >97% chromosomal ANI to *Y. pseudotuberculosis* |
| *C. botulinum* | BoNT cluster | Species is polyphyletic; cluster also in *C. butyricum* |
| STEC | stx1 / stx2 / eae | *E. coli* presence is commensal and uninformative |
| *V. cholerae* | ctxAB | Non-toxigenic environmental strains are widespread |
| *C. diphtheriae* | tox | Phage-borne, not chromosomal |
| *B. pertussis* | ptxA, not IS481 | IS481 is multi-copy and shared with *B. holmesii* |
| MTBC | RD9 present / RD1 absent | Members are ~99.9% identical; RD1 absence means BCG vaccine |
| *T. brucei* | TgsGP / SRA | *T. b. brucei* is not human-infective |
| *A. fumigatus* | cyp51A alleles | Presence is meaningless; resistance is the finding |

Three complementary evidence tiers run per sample:

- **dna** — read alignment to marker nucleotide sequence. Precise, blind to a divergent
  allele.
- **protein** — translated search (MMseqs2, high sensitivity). Sees through synonymous
  divergence.
- **hmm** — profile search. The remote-homology tier, and the reason a *novel* BoNT or Stx
  subtype, or a divergent filovirus gene, is still caught. This is what separates detecting
  a known strain from detecting a new agent.

A marker is PRESENT if **any** tier passes. They are complementary, not a consensus vote —
requiring agreement would discard exactly the sensitivity the protein and HMM tiers exist
to provide.

### Three-valued logic, and why it matters

Rules are boolean expressions over marker IDs, evaluated under **Kleene three-valued
logic**, because a marker has three states and not two:

```
PRESENT   searched for, found above threshold
ABSENT    searched for, confidently not found
UNKNOWN   never actually searched
```

"Absent" is evidence. "Untested" is not. Collapsing them lets a gap in the marker database
read as a negative result, so the pipeline would happily conclude *"RD1 absent, therefore
BCG"* for a sample it never tested, or *"no pXO2, therefore not anthrax"* for a marker it
never had.

The consequence that does the work: `UNKNOWN and FALSE` is FALSE — one confident failure
decides — but `UNKNOWN and TRUE` is UNKNOWN, because you cannot claim success on untested
evidence. **A rule evaluating to UNKNOWN never promotes a finding.** It holds it at
`supported` and names the untested markers.

`scripts/build_marker_db.py` writes a manifest recording, per marker, whether usable
sequence was actually obtained. That manifest — not the absence of alignment hits — is what
makes ABSENT calls legitimate.

Rules are parsed by a recursive-descent parser, never `eval()`, so a rules file cannot
execute code. Operators: `and or not present(x) absent(x) any_of(...) all_of(...)
n_of(k, ...)`.

```yaml
"Bacillus anthracis":
    confirm_requires: (banth_pxo1_pagA or banth_pxo1_lef or banth_pxo1_cya) and (banth_pxo2_capB or banth_pxo2_capA)
    veto_if: bcereus_chr_ces
```

**Identical** alignment evidence (breadth 0.38, occupancy 0.60, 900 reads, beating
neighbours 64:1) produces four different outcomes depending only on markers:

| Marker state | Tier |
|---|---|
| pXO1 + pXO2 present | `confirmed_in_silico` |
| pXO1 only, pXO2 confidently absent | `supported` — required markers not found |
| pXO1 only, pXO2 **untested** | `supported` — marker evidence incomplete (database gap) |
| pXO1 + pXO2 present, *B. cereus* cereulide also present | `marker_veto` |

And MTBC with rpoB present and RD1 absent → `vaccine_strain` (BCG), not disease.

Markers are **keyed by gene symbol, not accession**, for the same reason the watchlist is
keyed by name: accession numbers written from memory are a fabrication risk in a file that
gates a high-consequence call. `build_marker_db.py` resolves each row to real sequence and
records the accession, source and retrieval date.

Coverage: 166 markers, 59 rule sets over 90 watchlist taxa. **Every `high_consequence`
*species* has a rule.** The one `high_consequence` entry without one is `Filoviridae`, a
deliberate group-rank row capped at `group_lead`, so a marker rule would have nothing to
gate. *B. mallei* deliberately has no `confirm_requires` — it is a clonal deletion
derivative defined by ABSENT regions, which presence markers cannot express, so it caps at
`supported` rather than producing an unsupportable call.

An unevaluable **vaccine or veto** check now caps the tier at `supported`. This reduces
confirmations deliberately: "supported, vaccine status untestable" costs a follow-up,
"confirmed wild-type poliovirus" costs an outbreak response.

## Reference masking: the specificity lever

Two independent 2025–2026 reviews reached the same conclusion: **curating and masking the
reference database buys more specificity than anything available downstream, including
assembly**, at almost no sensitivity cost — 40–94% fewer false detections with 100%
true-positive preservation (*J Clin Microbiol* 2026).

The concrete failure mode, from Wright et al. 2026: an *H. pylori* call had **zero** reads
mapping to its reference once checked, traced to contaminating **phage sequence inside the
reference genome**. An unmasked reference does not merely fail to help — it manufactures
false positives no downstream threshold can catch, because the reads genuinely do align.

Four categories are identified:

| Source | Why |
|---|---|
| `lowcomplexity` | Homopolymers, tandem repeats, low-complexity windows — the classic spurious-mapping substrate |
| `vector` | NCBI UniVec_Core. Cloning-vector and adapter sequence left in assemblies is pure artifact |
| `host` | The pack's own host genome. A host-derived region in a pathogen reference collects the host reads that survived depletion |
| `crossgroup` | Regions also present near-identically in a **phylogenetically distant** reference in the same pack. Same-group hits are deliberately *not* masked — shared sequence between close relatives is real biology, and resolving it is what competitive mapping is for |

### Soft-mask, not hard-mask

Masking records **BED intervals with a reason**; it does not overwrite sequence with Ns.
Downstream, `coverage_stats.py --mask-bed` excludes masked bases from the breadth and
occupancy **denominators**, excludes reads lying entirely inside masked intervals from the
read count — and still reports them as `reads_in_masked_regions`, so the situation is
visible.

Hard-masking would make the Wright case *invisible*: the reads would stop aligning and the
taxon would quietly drop out with no record of why. A finding whose support is dominated by
masked-region reads is capped at `supported` with a caveat naming the cause.

If masking would exceed `--max-masked-fraction` (default 20%), the genome is flagged
`OVER_MASKED` rather than applied silently. A reference that is 60% masked is not curated,
it is misconfigured.

### What it actually found

Over the built human pack — 375 genomes, 5.76 Gbp:

| Source | Masked | Spread |
|---|---|---|
| low complexity | 19,864,369 bp | the bulk, as expected |
| **vector** | **104,184 bp** | **45 contigs** |
| crossgroup | 81,196 bp | 50 contigs (bacteria 7, fungus 12, parasite 31, virus 0) |
| host | 2,737 bp | 5 contigs |
| **total** | **20,040,354 bp (0.348%)** | no genome tripped `OVER_MASKED` |

Per-source rows are **upper bounds, not a partition** — an interval masked for two reasons
counts under both.

Read that honestly: 0.348% is nowhere near the headline 40–94% reduction, and it should not
be. 297 of 375 references here are RefSeq reference-grade complete assemblies, already well
curated — the *Y. pestis* cross-group run correctly found **0 bp**. What matters is the
104.4 kb of cloning-vector sequence in **45 reference contigs** and 81.9 kb of
cross-kingdom sequence in 50 more: exactly the mechanism behind the Wright false positive,
now excluded from every breadth denominator with the reads reported rather than silently
counted. Expect the payoff to concentrate in the 40 scaffold- and contig-level assemblies
relaxed to for several fungi and parasites — the ragged tail is where an unnoticed false
positive would come from. Virus cross-group masking found 0 bp, which is the right answer.

## Domain packs

One database for all three domains would be unusable. The fix is not a bigger index, it is
**scope**: one shared general-purpose Kraken2 index answers "what is in here", and
per-domain **packs** supply what is actually domain-specific — the host genome for
depletion, the watchlist, curated reference genomes, near-neighbour sets, and per-taxon
thresholds.

A pack is **data, not code**. Adding a domain or retargeting one means editing
`packs/<name>/pack.yaml` and its watchlist. The pipeline never changes.

```
packs/human/
  pack.yaml                  host, DB paths, thresholds, disclaimer
  watchlist.tsv              curated, keyed by NAME
  watchlist.resolved.tsv     generated: names -> taxids (+ descendants)
  refs.lock.tsv              every reference assembly, pinned by accession and checksum
```

### Watchlists are keyed by name, not taxid

NCBI merges, splits and retires taxonomy nodes constantly. A hardcoded taxid rots silently,
and a rotted taxid means an organism stops being screened for with **no error** — the worst
available failure mode here.

So watchlists carry names, and `bin/resolve_watchlist.py` reconciles them against a pinned
local taxdump at build time. Anything that does not resolve to exactly one node is a **hard
error** with candidates printed; the pipeline refuses to run on an unresolved watchlist.
That paid off on the first run: it caught 18 virus entries whose NCBI rank had changed, 11
taxa ICTV has since renamed, and 4 malformed near-neighbour entries.

It also expands each watchlist taxon to all its **descendants** — 90 seed taxa expand to
~134,500 taxonomy nodes. Not decoration: Kraken2 assigns many reads to *internal* nodes, so
matching only the species taxid misses most of the signal.

### The seed watchlist is a starting point

`packs/human/watchlist.tsv` covers organisms appearing on published notifiable-disease and
high-consequence lists, with per-taxon thresholds and caveats. **It must be reconciled
against the authoritative list for your jurisdiction before operational use** — for Canada
the PHAC ePATHogen risk-group database and the national notifiable disease list; WOAH
listed diseases for animal; EPPO A1/A2 for plant.

Every watchlist row carries a `notes` field that travels into the report. These are the
traps that make naive metagenomic pathogen calling dangerous:

- *B. anthracis* chromosome is near-identical to the B. cereus group — a species call needs
  plasmid markers, not chromosomal breadth.
- *Brucella* species are >99% ANI to each other; species-level calls from metagenomes are
  usually not defensible.
- *Entamoeba dispar* is morphologically identical to *E. histolytica* and non-pathogenic —
  the single biggest false positive in that genus.
- *Burkholderia thailandensis* is common soil flora and the main false-positive route to
  *B. pseudomallei*.
- Reston ebolavirus is not pathogenic to humans, so species resolution within *Ebolavirus*
  is essential rather than optional.
- Vaccine and attenuated strains (vaccinia, BCG, Sabin, 17D, RB51, TC-83) are routine
  confounders needing explicit discrimination from wild type.

## Comparing tools rather than arguing about them

Which classifier is best for this job is an empirical question, so the harness answers it.
`scripts/run_bakeoff.sh` runs the pipeline once per configuration over identical simulated
data and compares accuracy **and cost** — CPU-minutes and peak RSS, harvested from each
run's Nextflow trace. Variants are built as **minimal pairs**: a variant that changes two
things at once cannot attribute the difference to either. See
[`docs/benchmarking.md`](https://github.com/duceppemo/bioThreatFinder/blob/main/docs/benchmarking.md).

## References

The pipeline runs: [Nextflow](https://doi.org/10.1038/nbt.3820),
[fastp](https://doi.org/10.1093/bioinformatics/bty560),
[chopper](https://github.com/wdecoster/chopper),
[minimap2](https://doi.org/10.1093/bioinformatics/bty191),
[samtools](https://doi.org/10.1093/gigascience/giab008),
[Kraken2](https://doi.org/10.1186/s13059-019-1891-0),
[Bracken](https://doi.org/10.7717/peerj-cs.104),
[sylph](https://doi.org/10.1038/s41587-024-02412-y),
[MMseqs2](https://doi.org/10.1038/nbt.3988),
[HMMER](https://doi.org/10.1371/journal.pcbi.1002195),
[seqkit](https://doi.org/10.1371/journal.pone.0163962), and optionally
[Metabuli](https://doi.org/10.1038/s41592-024-02273-y) and
[KrakenUniq](https://doi.org/10.1186/s13059-018-1568-1).

Design evidence cited above: Oskolkov 2026 (*Front Microbiol*,
doi:10.3389/fmicb.2026.1603339) on Kraken2 versus KrakenUniq; Wright et al. 2026 on
reference contamination; *J Clin Microbiol* 2026 on masking as a specificity lever;
[Perseus](https://doi.org/10.1093/bioinformatics/btag687) on confident misassignment of
sequence whose family is absent from the database. Every choice, its evidence tag, and what
would reverse it is in
[`docs/decisions.md`](https://github.com/duceppemo/bioThreatFinder/blob/main/docs/decisions.md).
