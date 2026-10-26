<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/logo/bioThreatFinder-logo-dark.svg">
    <img alt="bioThreatFinder" src="assets/logo/bioThreatFinder-logo-light.svg" width="440">
  </picture>
</h1>

<p align="center">
  <a href="https://github.com/duceppemo/bioThreatFinder/actions/workflows/tests.yml"><img src="https://github.com/duceppemo/bioThreatFinder/actions/workflows/tests.yml/badge.svg" alt="tests"></a>
  <a href="https://github.com/duceppemo/bioThreatFinder/actions/workflows/containers.yml"><img src="https://github.com/duceppemo/bioThreatFinder/actions/workflows/containers.yml/badge.svg" alt="containers"></a>
  <a href="https://github.com/duceppemo/bioThreatFinder/releases/latest"><img src="https://img.shields.io/github/v/release/duceppemo/bioThreatFinder?label=release&cacheSeconds=3600" alt="Latest release"></a>
  <img src="https://img.shields.io/badge/nextflow-%E2%89%A524.04-23aa62" alt="Nextflow >=24.04">
  <a href="LICENSE"><img src="https://img.shields.io/github/license/duceppemo/bioThreatFinder" alt="License: MIT"></a>
  <a href="https://github.com/duceppemo/bioThreatFinder/wiki"><img src="https://img.shields.io/badge/docs-wiki-informational" alt="Documentation"></a>
  <a href="https://doi.org/10.5281/zenodo.23607033"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.23607033.svg" alt="Software DOI"></a>
  <a href="https://doi.org/10.5281/zenodo.23602998"><img src="https://img.shields.io/badge/databases-10.5281%2Fzenodo.23602998-blue" alt="Databases DOI"></a>
</p>

Metagenomic screening and in-silico confirmation of high-consequence pathogens, across
human, animal and plant health. A permissive stage-1 screen across several classifier arms
feeds a strict, per-candidate stage 2 that maps reads against each target **together with
its near neighbours**, so reads that fit a harmless relative better are attributed to the
relative rather than to the agent. A marker-gene tier gates the top evidence tier wherever
the chromosome is not the actionable signal — pXO1/pXO2 for anthrax, not chromosomal
breadth.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/workflow/bioThreatFinder-workflow-dark.svg">
  <img alt="bioThreatFinder workflow: a permissive stage-1 OR across classifier arms feeds a strict, per-candidate stage 2 that maps reads against the target together with its near neighbours, so reads preferring a harmless relative are attributed to the relative." src="assets/workflow/bioThreatFinder-workflow-light.svg" width="100%">
</picture>

> **Research software, not a validated diagnostic assay.** Sequence detection does not
> establish a viable organism, an active infection, or exposure. No threshold here is a
> measured operating point, and **nothing this produces should be used to decide that a
> sample is safe.** Every finding needs orthogonal laboratory confirmation and
> interpretation by a qualified laboratorian. Read [`DISCLAIMER.md`](DISCLAIMER.md) first.

## Quick start

```bash
# Tools: three conda environments, or containers (-profile docker|podman|singularity|apptainer)
conda env create -f envs/btf-refs.yml      # and btf-kraken2.yml, btf-sylph.yml

export BTF_DB=/fast/disk/bioThreatFinder   # must end in bioThreatFinder
export BTF_DB_ROOT=/shared/read-only       # pre-existing shared resources

scripts/install_databases.sh --route minimal     # ~62 GB, sylph stage 1
nextflow run . -profile conda --input samplesheet.csv --pack human --outdir results
```

The installer downloads only the two artifacts that **cannot be re-created** — a pinned
NCBI taxonomy snapshot and the curated marker database, 166 MB in total
([doi:10.5281/zenodo.23602998](https://doi.org/10.5281/zenodo.23602998)) — verifying each
sha256 before unpacking, then builds the rest. The 43 GB reference pack is built from
[`packs/human/refs.lock.tsv`](packs/human/refs.lock.tsv), which pins every assembly
accession **and** the checksum of its decompressed sequence, so a clone builds the pack the
measurements were made on and can prove it. `--route standard` adds the 115 GB Kraken2
index — which needs **~700 GB free while it builds**, several times its finished size;
`--route full` adds Metabuli.

```csv
sample,platform,fastq_1,fastq_2,role
SPEC001,illumina,/data/S1_R1.fastq.gz,/data/S1_R2.fastq.gz,sample
SPEC002,ont,/data/S2.fastq.gz,,sample
NTC01,illumina,/data/NTC_R1.fastq.gz,/data/NTC_R2.fastq.gz,control
```

**Include a negative control.** Without one the pipeline cannot subtract reagent,
index-hopping or carry-over contamination, and it stamps a run-level warning on every
finding. Open `results/report/biothreatfinder_report.html` for the run, or
`results/score_evidence/findings.tsv` for the machine-readable result — every finding
carries the gate that decided it, never a bare verdict.

```bash
python3 tests/test_pipeline.py     # 517 checks with databases; 451 on a bare clone, ~30 s
```

## Documentation

Everything else is in the [wiki](https://github.com/duceppemo/bioThreatFinder/wiki), whose
sources are maintained in [`docs/wiki`](docs/wiki):

- [Installation](https://github.com/duceppemo/bioThreatFinder/wiki/Installation) — the three environments and why three, the two database roots, containers
- [Usage](https://github.com/duceppemo/bioThreatFinder/wiki/Usage) — samplesheet, options, profiles, building each database by hand
- [Methods](https://github.com/duceppemo/bioThreatFinder/wiki/Methods) — why two stages, the evidence tiers, the classifier arms, the marker tier's three-valued logic, reference masking
- [Databases](https://github.com/duceppemo/bioThreatFinder/wiki/Databases) — what each artifact is, what is published and what is built, and the checks that answer whether one still covers its inputs
- [Outputs](https://github.com/duceppemo/bioThreatFinder/wiki/Outputs) — files, columns, the stage-1 coverage table, Pavian
- [Validation](https://github.com/duceppemo/bioThreatFinder/wiki/Validation) — what has been measured, and the specific gaps
- [FAQ](https://github.com/duceppemo/bioThreatFinder/wiki/FAQ) — empty results, exit code 3, and the other things that look like bugs
- [Development](https://github.com/duceppemo/bioThreatFinder/wiki/Development) — tests, CI, releases, and the two conventions that are not optional

The engineering record stays in the repository, because it is written for whoever maintains
this rather than for whoever runs it:
[`docs/decisions.md`](docs/decisions.md) (why each tool was chosen, the evidence, and what
would reverse it — including the decisions that were reversed),
[`docs/operations.md`](docs/operations.md),
[`docs/verification.md`](docs/verification.md),
[`docs/benchmarking.md`](docs/benchmarking.md),
[`docs/validation.md`](docs/validation.md),
[`docs/samplesheet.md`](docs/samplesheet.md),
[`docs/roadmap.md`](docs/roadmap.md).

## Citation

If bioThreatFinder is useful in your work, please cite the software (see
[`CITATION.cff`](CITATION.cff)) and the databases:

> Duceppe, M.-O. (2026). *bioThreatFinder: two-stage metagenomic screening and in-silico
> confirmation of high-consequence pathogens* (v0.4.1). Zenodo.
> https://doi.org/10.5281/zenodo.23607033

> Duceppe, M.-O. (2026). *bioThreatFinder databases: a pinned NCBI taxonomy snapshot
> (2026-09-29) and a curated marker-gene database for high-consequence pathogen
> screening*. Zenodo. https://doi.org/10.5281/zenodo.23602998

and the programs it runs: Nextflow, fastp or chopper, minimap2, samtools, Kraken2, Bracken,
sylph, MMseqs2 and HMMER, with Metabuli and KrakenUniq when enabled
([references](https://github.com/duceppemo/bioThreatFinder/wiki/Methods#references)).

## License

[MIT](LICENSE). The published databases are derived data under
[CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/); see
[`DISCLAIMER.md`](DISCLAIMER.md) for what they are and are not.
