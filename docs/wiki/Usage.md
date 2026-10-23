# Usage

```bash
nextflow run . -profile conda,workstation \
    --input samplesheet.csv --pack human --outdir results
```

## The samplesheet

```csv
sample,platform,fastq_1,fastq_2,role
SPEC001,illumina,/data/S1_R1.fastq.gz,/data/S1_R2.fastq.gz,sample
SPEC002,ont,/data/S2.fastq.gz,,sample
NTC01,illumina,/data/NTC_R1.fastq.gz,/data/NTC_R2.fastq.gz,control
```

`platform` is `illumina` or `ont`; `role` is `sample` or `control`. An empty samplesheet is
a hard error, not a green run that analysed nothing.

**Include a negative control.** Without one the pipeline cannot subtract reagent,
index-hopping or carry-over contamination, and it stamps a run-level warning on every
finding in the report. Kit contaminants and barcode bleed are the norm in low-biomass
metagenomics, not the exception.

What `role=control` actually does is worth reading before trusting a negative: control
subtraction is **run-global and keyed on taxid alone**, which is deliberate — a matrix
blank is supposed to reach the whole run — and it means one contaminated control can
suppress another scenario's genuine positive. The column-by-column account, including how
to read the audit table, is in
[`docs/samplesheet.md`](https://github.com/duceppemo/bioThreatFinder/blob/main/docs/samplesheet.md).

## Profiles

| Profile | What it does |
|---|---|
| `conda` / `mamba` | resolve each process's conda directive |
| `docker` / `podman` / `singularity` / `apptainer` | one image per process ([Installation](Installation)) |
| `workstation` | an example sizing for a single large machine; set your own |
| `test` | tiny inputs, heavy steps skipped — **does not stand alone**, see below |
| `debug` | keep work directories, dump hashes, echo the hostname |

`-profile test` exercises the workflow; it is not a self-contained test and cannot be made
into one by adding a samplesheet, because stage 1 needs a classifier database and stage 2
the reference pack — 60 GB that cannot live in a git repository. It says so at launch. For
a test that needs nothing at all, run `python3 tests/test_pipeline.py`.

## Options worth knowing

**The opt-in arms take an explicit value.** `--skip_metabuli false` enables Metabuli;
omitting the flag leaves it off. On the Nextflow CLI a bare `--flag false` arrives as the
string `"false"`, which is truthy in Groovy, so these spellings were silent no-ops until
`lib/BtfBool.groovy` began coercing them. If a flag seems to do nothing, that was why.

| Flag | Default | Effect |
|---|---|---|
| `--skip_kraken2` / `--skip_sylph` | false | the two arms that run by default |
| `--skip_metabuli` | **true** | DNA + amino-acid arm; has changed no finding on any scenario measured, kept for the case that cannot be simulated |
| `--skip_krakenuniq` | **true** | distinct-k-mer second opinion; needs a separate index |
| `--skip_markers` | false | the marker tier that gates `confirmed_in_silico` |
| `--skip_host_depletion` | false | host removal before everything else |
| `--publish_reads` | false | trimmed FASTQs are **not** published by default: they still contain host sequence |
| `--pavian_html` | false | render Pavian's static report; host environment only |
| `--contest_attribution` | false | enables `indistinguishable_complex`; parked pending recalibration |
| `--high_unclassified` | 0.40 | flag threshold for the unclassified fraction; uncalibrated, set from your own negatives |

Database paths come from `pack.yaml`; `--kraken2_db`, `--sylph_db`, `--metabuli_db` and
`--krakenuniq_db` are explicit overrides that default to null.

## Building the databases by hand

`scripts/install_databases.sh` runs these in order. Each is resumable and each refuses
rather than guessing when a root is unset.

```bash
# 1. Pin a taxonomy snapshot, then resolve the watchlist against it.
#    resolve_watchlist.py refuses a snapshot older than 90 days with
#    --require-fresh-taxdump, because a stale one silently loses clade coverage.
TAX="$BTF_DB/taxdump/$(date -u +%Y-%m-%d)"
mkdir -p "$TAX" && curl -sSo - \
  https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xz -C "$TAX"
python3 bin/resolve_watchlist.py \
    --watchlist packs/human/watchlist.tsv \
    --taxdump "$TAX" --require-fresh-taxdump \
    --out packs/human/watchlist.resolved.tsv \
    --descendants-out packs/human/watchlist.descendants.tsv

# 2. The databases. --taxdump must be the snapshot the watchlist was resolved against,
#    and the scripts refuse to run if it is not.
conda activate btf-refs
bash    scripts/fetch_host_genome.sh human --threads 24
python3 scripts/build_pack_refs.py  --pack human --from-lock packs/human/refs.lock.tsv \
                                    --taxdump "$TAX" --threads 24
python3 scripts/mask_pack_refs.py   --pack human --threads 24
python3 scripts/build_marker_db.py  --pack human --taxdump "$TAX" --threads 24

conda activate btf-sylph
bash    scripts/build_sylph_db.sh   --pack human --threads 24

conda activate btf-kraken2
bash    scripts/build_kraken2_db.sh --taxdump "$TAX" --threads 40
```

**`--from-lock` is not optional if you want the pack the measurements were made on.**
Without it each taxon resolves to whatever NCBI calls latest that day: the rule is
deterministic, the catalogue is not, so two people build two different packs from the same
watchlist with no error on either side. See [Databases](Databases).

**Order matters.** The sylph sketch records the genome set it covers, so one built before
the references silently omits them; masking needs the host genome for its host pass.

## Resource sizing

`conf/base.config` sets per-process requests and clamps them with `resourceLimits`. The
`workstation` profile is an example for one large machine (60 CPUs, 420 GB) — set your own,
or leave `resourceLimits` to do the clamping. On a cluster, add an executor in a site
config passed with `-c`; Nextflow handles SLURM and PBS natively, and that is what makes
the translated-search tier affordable.

One timing fact worth planning around: `MARKER_SEARCH_AA` took **7h 20m** on a single ONT
sample in the first real long-read run, so the process time limit is 24 h rather than the
4 h it started at.
