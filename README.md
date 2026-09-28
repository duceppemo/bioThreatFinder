# bioThreatFinder

Metagenomic screening and in-silico confirmation of high-consequence pathogens,
across human, animal and plant health.

> **This tool produces analytical evidence, not diagnoses.** Sequence detection
> does not establish viable organism, active infection, or exposure. Every
> finding requires orthogonal laboratory confirmation and interpretation by a
> qualified laboratorian before any clinical or public-health action.

---

## Why two stages

Stage-1 metagenomic classifiers (Kraken2, sylph, Centrifuge) are built for
*community profiling*. They answer "what is the composition of this sample",
and they are good at it. They are not built to answer "is organism X present",
which is a fundamentally different question with a much harsher false-positive
cost.

A handful of k-mers shared with a common environmental relative is enough to
place a high-consequence taxon in a Kraken2 report. Run a classifier over a
soil or sputum metagenome and you will see plenty of alarming names in the
output. Essentially all of them are noise.

So the pipeline splits the question:

| Stage | Question | Tooling | Posture |
|---|---|---|---|
| 1 | What *might* be here? | Kraken2 (+Bracken), sylph | Deliberately permissive — missing a lead is worse than carrying a false one forward |
| 2 | Does the alignment evidence hold up? | Competitive minimap2 + coverage geometry | Strict — this is what decides |

Stage 2 is where the tool earns its keep. For every watchlist candidate it maps
reads against a bundle containing the **target plus its near neighbours**, then
asks four things that k-mer counting cannot:

- **Breadth** — what fraction of the reference is covered at all? A real
  organism at low abundance spreads thin across the whole genome. An artifact
  stacks on one conserved gene.
- **Window occupancy** — breadth's robust sibling. At very low depth, breadth is
  small either way, but 40 reads across 40 windows is an organism and 40 reads
  in 2 windows is a repeat.
- **Competitive attribution** — do these reads actually prefer the target, or do
  they prefer its harmless relative? `neighbour_ratio` below ~1 means the answer
  is the relative.
- **Control subtraction** — is this above what the negative control shows?

## Evidence tiers

Every finding is assigned a tier and carries a `reasons` string explaining which
gate decided it. Nothing is a bare verdict.

| Tier | Meaning |
|---|---|
| `confirmed_in_silico` | Met every threshold. **Still requires wet-lab confirmation.** |
| `supported` | Real alignment evidence, below the confirmation bar |
| `no_reference` | On the watchlist, raised by the screen, but no reference genome available — **a database gap, not a negative result** |
| `suppressed_control` | Comparable level in a negative control; contamination until proven otherwise |
| `attributed_to_neighbour` | Reads align better to a relative |
| `control_background` | What the control itself contains (QC, not a finding) |
| `screened_only` | Stage-1 lead that failed stage-2 — the expected fate of classifier noise |

## Domain packs

You were right that one database for all three domains would be unusable. The
fix is not a bigger index, it is **scope**:

- One **shared** general-purpose Kraken2 index answers "what is in here".
- Per-domain **packs** supply what is actually domain-specific: the host genome
  for depletion, the watchlist, curated reference genomes, near-neighbour sets,
  and per-taxon thresholds.

A pack is **data, not code**. Adding a domain or retargeting one means editing
`packs/<name>/pack.yaml` and its watchlist. The pipeline never changes.

```
packs/human/
  pack.yaml                  host, DB paths, thresholds, disclaimer
  watchlist.tsv              curated, keyed by NAME
  watchlist.resolved.tsv     generated: names -> taxids (+ descendants)
```

### Watchlists are keyed by name, not taxid

This is deliberate and it matters. NCBI merges, splits and retires taxonomy
nodes constantly. A hardcoded taxid rots silently, and a rotted taxid means an
organism stops being screened for with **no error** — the worst possible failure
mode here.

So watchlists carry names, and `scripts/resolve_watchlist.py` reconciles them
against a pinned local taxdump at build time. Anything that does not resolve to
exactly one node is a **hard error** with candidates printed. The pipeline
refuses to run on an unresolved watchlist.

This paid off on the first run: it caught 18 virus entries whose NCBI rank had
changed, 11 taxa that ICTV has since renamed (resolved via equivalent names),
and 4 malformed near-neighbour entries.

It also expands each watchlist taxon to all its **descendants** — 89 seed taxa
expand to ~128,000 taxonomy nodes. This is not decoration: Kraken2 assigns many
reads to *internal* nodes, so matching only the species taxid misses most of the
signal.

### The seed watchlist is a starting point

`packs/human/watchlist.tsv` covers organisms that appear on published
notifiable-disease and high-consequence-pathogen lists, with per-taxon
thresholds and detection caveats. **It must be reconciled against the
authoritative list for your jurisdiction before operational use** — for Canada
that means the PHAC ePATHogen risk-group database and the national notifiable
disease list; WOAH listed diseases for the animal pack; EPPO A1/A2 for plant.
Those lists change, and provenance matters for anything you would act on, so
they belong in a documented fetch step rather than baked into a file by whoever
wrote the tool.

### Taxon caveats are first-class

Every watchlist row carries a `notes` field that travels all the way into the
report. These are the traps that make naive metagenomic pathogen calling
dangerous:

- *Bacillus anthracis* chromosome is near-identical to the B. cereus group — a
  species call needs plasmid markers, not chromosomal breadth.
- *Brucella* species are >99% ANI to each other; species-level calls from
  metagenomes are usually not defensible.
- *Entamoeba dispar* is morphologically identical to *E. histolytica* and
  non-pathogenic — the single biggest false positive in that genus.
- *Burkholderia thailandensis* is common soil flora and the main false-positive
  route to *B. pseudomallei*.
- Reston ebolavirus is not pathogenic to humans, so species resolution within
  *Ebolavirus* is essential rather than optional.
- Vaccine and attenuated strains (vaccinia, BCG, Sabin, 17D, RB51, TC-83) are
  routine confounders and need explicit discrimination from wild type.

## Quick start

```bash
# 1. Resolve the watchlist (fast, and validates your taxonomy)
python3 scripts/resolve_watchlist.py \
    --watchlist packs/human/watchlist.tsv \
    --taxdump /media/36tb/db/taxdump \
    --out packs/human/watchlist.resolved.tsv \
    --descendants-out packs/human/watchlist.descendants.tsv

# 2. Build databases (hours; see scripts/ for space requirements)
bash   scripts/fetch_host_genome.sh human
python3 scripts/build_pack_refs.py --pack human --threads 32
bash   scripts/build_kraken2_db.sh --threads 48

# 3. Run
nextflow run . -profile conda,workstation \
    --input samplesheet.csv --pack human --outdir results
```

### Samplesheet

```csv
sample,platform,fastq_1,fastq_2,role
SPEC001,illumina,/data/S1_R1.fastq.gz,/data/S1_R2.fastq.gz,sample
SPEC002,ont,/data/S2.fastq.gz,,sample
NTC01,illumina,/data/NTC_R1.fastq.gz,/data/NTC_R2.fastq.gz,control
```

`platform` is `illumina` or `ont`. `role` is `sample` or `control`.

**Include a negative control.** Without one the pipeline cannot subtract
reagent, index-hopping or carry-over contamination, and it will stamp a
run-level warning on every finding in the report. Kit contaminants and barcode
bleed are the norm in low-biomass metagenomics, not the exception.

## Outputs

```
results/
  watchlist/        resolved watchlist + resolution log
  screen/           per-sample stage-1 candidates
  confirm/          per-(sample,taxon) coverage metrics
  score_evidence/   findings.tsv  <- the machine-readable result
  report/           biothreatfinder_report.html  <- the human-readable one
  pipeline_info/    timeline, trace, DAG
```

## Requirements

- **Nextflow >= 23.10** — the version currently on this box is 21.10.4 and
  needs upgrading
- conda/mamba (all tool envs are declared per process)
- ~200 GB RAM for a PlusPF-scale Kraken2 index; everything else is modest
- Analysis scripts in `bin/` are **stdlib-only Python 3** plus `samtools`, so
  they run in whatever env the pipeline hands them

## Status

v0.1.0 — scaffold complete and internally validated:

- Watchlist resolution tested against the live taxdump: 89/89 resolved
- Evidence tiering validated against synthetic fixtures covering every decision
  path (true positive, conserved-repeat artifact, neighbour attribution, control
  suppression, database gap, weak lead)
- Nextflow config and workflow graph parse

Not yet run on real sequencing data, and no databases built yet.

### Known gaps

- Animal and plant packs are directory stubs — no watchlist or host genomes
- Marker-gene and plasmid logic is specified in the watchlist notes but not yet
  implemented, so calls for *B. anthracis*, *C. botulinum*, STEC and toxigenic
  *V. cholerae* are not yet defensible
- Influenza A segment-aware subtyping is out of scope and needs a dedicated
  subtyper
- No AMR characterisation (CARD/ResFinder/DeepARG are already on this box and
  would slot in as a stage-3 module)
- Read-level output for candidate taxa is not extracted, so an analyst cannot
  eyeball the alignments yet
