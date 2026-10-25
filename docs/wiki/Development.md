# Development

Start with
[`CONTRIBUTING.md`](https://github.com/duceppemo/bioThreatFinder/blob/main/CONTRIBUTING.md),
which states the two conventions that are not optional. This page is the mechanics.

## The suite

```bash
python3 tests/test_pipeline.py                 # 483 checks, offline, no databases, ~30 s
python3 tests/test_pipeline.py --list          # the groups
python3 tests/test_pipeline.py -k markerdb     # one group — SUBSTRING match, not an expression
```

`-k` is a substring match. `-k "idxcover or selftests"` matches nothing, which has wasted
time before.

It is a self-contained runner, not pytest: fixture-driven, no network, no databases. Cases
that compare against a built database skip with a message naming what to set, so a bare
clone is a valid place to run it — 422 pass and 14 skip there. With `$BTF_DB` set and
`samtools` on `PATH`, 483 pass.

Per-script `--self-test` suites run as single cases from here, so `selftest:build_marker_db`
alone covers 60+ assertions.

The suite asserts **designed behaviour, not accuracy.** A failure means the logic changed;
it does not mean the pipeline got less accurate, because nothing here measures accuracy.

## What a new check has to do

**Mutation-test it.** Break the thing it is supposed to catch and watch it fail. Three
checks in this repository passed a mutation they should have caught and were rewritten —
one compared counts that were never meant to agree, one anchored a regex on `=` and missed
every default reached through a conditional, one exempted comments and so exempted exactly
where the problem was.

**And make sure it cannot fire on a correct installation.** Four guards were caught doing
that and fixed before shipping: one called the correct shipped database stale on 116 of 166
markers, one reported a replaced assembly for byte-identical sequence fetched by a
different route, one compared a directory against a literal Groovy expression and refused
every build on every machine, and one smoke-tested images through a login shell that drops
the environment's `PATH`. A guard that cries wolf gets switched off, which is worse than
not having it.

**Distinguish *no* from *could not tell*.** The database checks exit `0` for yes, `3` for no
with the reason named, and `1` for undetermined. Conflating the last two reports a missing
tool as a stale database.

## Generated files

Change the source, then regenerate — the suite fails if they drift:

| Generated | From | By |
|---|---|---|
| `envs/process/*.yml`, `conf/containers.config` | the `conda` directives in `modules/local/*.nf` | `scripts/sync_container_specs.py` |
| `packs/*/refs.lock.tsv` | a built reference pack | `build_pack_refs.py --write-lock ... --dry-run` |
| `packs/*/watchlist.resolved.tsv` | `watchlist.tsv` + the pinned taxdump | `bin/resolve_watchlist.py` |
| `assets/databases.manifest.tsv` | the packaged archives | `scripts/package_databases.sh` |
| the workflow SVG | `main.nf` | `scripts/make_workflow_svg.py` |

A hand-maintained container map is the one that rots: a pin moves, nobody updates the map,
and the conda route and the container route run different versions of the same tool —
correct one way, wrong the other, visible only in results.

## No absolute paths, anywhere

Not in code, not in a comment, not in documentation. `tests/test_pipeline.py -k roots`
scans **every tracked file** via `git ls-files`; the only exemption is `docs/decisions.md`,
a dated log whose entries cite the paths their measurements were taken on. Database
locations come from `$BTF_DB` / `$BTF_DB_ROOT`, and every script refuses rather than
defaulting.

## CI

| Workflow | When | What |
|---|---|---|
| `tests.yml` | push, PR | the suite on Python 3.11 and 3.13 **with samtools** (without it the suite passes while registering 19 fewer cases); compiles `main.nf`; checks the unset-roots refusal and that every build script refuses with no `$BTF_DB`; ruff and shellcheck, advisory |
| `tests.yml` (`conda-specs`) | weekly, manual | solves every per-process spec and the three build environments. Network-bound and minutes long, so not on every commit |
| `containers.yml` | version tag, PR touching container files | builds all 14 images, **smoke-tests every pinned tool against the version it reports**, pushes to ghcr.io on a tag, and pulls each image **with no credentials** afterwards |
| `wiki.yml` | push touching `docs/wiki/**` | syncs these pages to the GitHub wiki |

Compiling `main.nf` is not ceremony: a Groovy closure in a local `def` reached a committed
file once, which Nextflow 26's strict syntax refuses to compile — and markers are enabled
by default, so that was every run. Reading it did not catch that.

The anonymous-pull job exists because GHCR publishes every new package **private** by
default with no API to change it, so a release can end with fourteen images built, tested
and pushed while `-profile docker` fails for everyone else.

## Releases

1. Update `version` in `nextflow.config`'s manifest and in `CITATION.cff` — the suite
   checks they agree, because a DOI that names a version the software does not report is
   worse than no version.
2. Regenerate anything derived from the version: `scripts/sync_container_specs.py` writes
   the image tags literally.
3. `git tag -a vX.Y.Z` and `gh release create`. The tag triggers the container build and
   push.
4. Flip any **new** GHCR package to public (UI only, once per package).
5. If the databases changed, `scripts/package_databases.sh` and a new Zenodo version, then
   the urls in `assets/databases.manifest.tsv`.

The Zenodo deposit's metadata comes from `.zenodo.json`, which **overrides** `CITATION.cff`
for that purpose — so the two must agree, and the suite checks the author, licence and the
link to the database deposit. Without `.zenodo.json`, Zenodo infers from `CITATION.cff`,
which works and produces a software record with **no link to the databases the pipeline
cannot run without**.

Two DOIs, and the distinction matters when citing: the **concept** DOI always resolves to
the latest version and is what belongs in a paper; the **version** DOI pins exact files and
is what `assets/databases.manifest.tsv` uses, because the installer verifies a recorded
sha256 and "latest" would one day stop matching it.

| Record | Concept | Version |
|---|---|---|
| software | `10.5281/zenodo.23607033` | one per release |
| databases | `10.5281/zenodo.23602998` | `10.5281/zenodo.23602999` |

## Decisions

Every non-obvious choice is a numbered entry in
[`docs/decisions.md`](https://github.com/duceppemo/bioThreatFinder/blob/main/docs/decisions.md),
with what was decided, the evidence, and what would reverse it. Evidence is tagged
*independent / tool-paper / read-from-source / measured-here / reasoned*, so a reader can
tell a measurement from a judgement.

Several entries record decisions that were later **reversed by measurement**, and they were
left in rather than tidied away — including the one that made KrakenUniq primary, and the
one that chose a `.mailmap` over rewriting history. If your change contradicts an entry, say
so in the new one and explain what changed. Please do not quietly delete a finding that
turned out to be inconvenient.

## Layout

```
main.nf                 the workflow, plus the launch-time guards
modules/local/*.nf      processes, each with its own conda pin
conf/                   base.config (resources), containers.config (generated), test.config
bin/                    analysis scripts: stdlib-only Python 3 plus samtools
scripts/                database builders, checks, benchmarks, packaging
packs/<name>/            a domain pack: watchlist, markers, rules, thresholds, refs lock
envs/                   the three build environments + per-process container specs
containers/Dockerfile   one image per process spec
tests/test_pipeline.py  the whole suite
docs/                   the engineering record; docs/wiki/ is this wiki
```
