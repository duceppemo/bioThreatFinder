# Installation

## Requirements

- **Nextflow ≥ 24.04.** `conf/base.config` uses `resourceLimits`, which arrived in 24.04,
  and the config and DSL were migrated to the strict syntax Nextflow 24+ enforces — so an
  older Nextflow passes a loose manifest check, starts, and then fails to *parse*, which
  reads as a broken pipeline rather than an old Nextflow. Developed and run on 26.04.6.
  Note that migration rejects `check_max()` helpers, top-level `def` in config, and
  top-level `workflow.onComplete`, so any pre-24 pipeline on the same machine will also
  stop parsing.
- **Tools:** either three conda environments (below) or containers. Both routes read the
  same pins.
- **~200 GB RAM** for a PlusPF-scale Kraken2 index. MMseqs2 sensitivity settings dominate
  CPU time. Without the Kraken2 arm (`--route minimal`), requirements are far lower.
- **Disk:** ~62 GB for the minimal route, ~323 GB with Kraken2, ~433 GB with Metabuli.
  Fast local disk matters: the k-mer indexes are read-heavy and load time dominates
  per-sample wall clock.
- Analysis scripts in `bin/` are **stdlib-only Python 3** plus `samtools`, so they run in
  whatever environment the pipeline hands them.

## Tools

### Conda — three environments, and the split is not tidiness

```bash
conda env create -f envs/btf-refs.yml      # datasets, minimap2, seqkit, samtools, mmseqs2, hmmer, mafft
conda env create -f envs/btf-kraken2.yml   # kraken2 2.17.2, bracken
conda env create -f envs/btf-sylph.yml     # sylph 1.0.0
```

**Why kraken2 is alone.** Installing it alongside mmseqs2/hmmer/mafft makes the solver
resolve **kraken2 2.0.8-beta (2019)**, which has no `--report-minimizer-data` — the flag
the whole stage-1 design rests on. An index built with it *works* and carries none of the
evidence the pipeline scores on. `build_kraken2_db.sh` now refuses anything below 2.1.0 and
verifies the flag is actually advertised. The environment files are hand-curated rather
than `conda env export` output for the same reason: an export of the development machine's
own `btf-refs` would ship that trap as the recommended install.

**Why sylph is alone and pinned.** A sketch must be readable by the version that queries
it, and `modules/local/detect.nf` pins `sylph=1.0.0`. The 0.6 → 1.0 jump was a major
version change, so a sketch built by another major version is a database the run cannot
open.

The per-process pins are separate from these three and live in `modules/local/*.nf`;
`-profile conda` resolves them per task. `scripts/check_conda_specs.sh` puts every distinct
spec through the solver — it needs the network and takes minutes, so it is not in the unit
suite and runs weekly in CI instead.

### Containers

```bash
nextflow run . -profile docker ...        # or podman, singularity, apptainer
```

Fourteen small images, one per distinct process spec, at
`ghcr.io/duceppemo/biothreatfinder-<spec>:<version>`. They are **generated from the same
`conda` directives** the conda route uses — `envs/process/*.yml` and
`conf/containers.config` are written by `scripts/sync_container_specs.py`, and
`--check` fails the suite if they drift from the modules. A hand-maintained container map
is the one that rots: a pin moves, nobody updates the map, and the two routes run different
versions of the same tool — correct one way, wrong the other, visible only in results.

**Not one big image**, because that reproduces the kraken2 solver trap above. Per-spec
images also build in parallel, and a pin change rebuilds one.

```bash
# build locally instead of pulling, then point the pipeline at the local tags
docker build -f containers/Dockerfile --build-arg SPEC=sylph -t biothreatfinder-sylph:0.4.0 .
nextflow run . -profile docker --container_prefix biothreatfinder- ...
```

**`--pavian_html` is refused under any container engine.** Pavian and two of its
dependencies install from GitHub, so no conda spec can declare them and no image contains
them; that process runs in the host environment. The refusal happens at launch rather than
inside an image with no R.

## The two database roots

```bash
export BTF_DB_ROOT=/shared/read-only              # or --db_root
export BTF_DB=/fast/disk/bioThreatFinder          # or --btf_db
```

| Root | Wants | Holds |
|---|---|---|
| `$BTF_DB_ROOT` | capacity; may be slow and shared | Pre-existing **read-only** resources a site already has: a taxdump, a sylph GTDB sketch, `nt` |
| `$BTF_DB` | **the fastest disk available**, ~350 GB for everything | Project-owned: pack references and markers, the sylph **pack** sketch, host genomes, and the Kraken2 / Metabuli / KrakenUniq indexes |

Both must be set. The pipeline refuses to start with either unset, naming them, and **no
script in this repository defaults to a path on any particular machine** — they used to
default to the development host's two volumes, so a clone built against paths it did not
have. A test scans every tracked file to keep it that way.

**`$BTF_DB` must end in `bioThreatFinder`.** The build scripts that take `--db-root` append
that segment themselves, so any other name puts the builders and the pipeline in two
different trees — and references under a base the pipeline does not read make every
watchlist taxon report `no_reference`, which reads as "checked and clear" rather than
"never built".

Two flags override, and they do **not** mean the same thing: `--btf-db` is `$BTF_DB`
itself, while `--db-root` is its **parent**. `--db-root` is not `$BTF_DB_ROOT`; passing the
shared root to it is the mistake it used to default to. Neither flag is needed once the
variables are exported.

## Databases

```bash
scripts/install_databases.sh --route minimal     # ~62 GB
scripts/install_databases.sh --route standard    # + Kraken2 PlusPF, ~323 GB
scripts/install_databases.sh --route full        # + Metabuli, ~433 GB
```

Resumable — rerun the same command after an interruption. See [Databases](Databases) for
what each route contains, what is downloaded versus built, and the read-only checks that
answer whether a database still covers its inputs.

**`minimal` is not small and not crippled.** Stage 2 maps against the reference pack and
sylph profiles it, so neither the pack nor the host genome is optional. What it leaves out
is the 261 GB Kraken2 index — and 19 of 90 watchlist taxa have no node in that index at
all, so sylph is already the only stage-1 arm that can see them.

## Checking an installation

```bash
python3 tests/test_pipeline.py              # 477 checks, offline, no databases, ~30 s
nextflow run . --help                       # must compile
```

On a bare clone with no databases and no samtools, 416 of those checks pass and 14 skip
with a message naming what to set — so the suite is a valid thing to run before anything
is built. Then the database checks in [Databases](Databases).
