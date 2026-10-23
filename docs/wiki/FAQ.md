# FAQ

Mostly: things that look like bugs and are not, and things that look fine and are not.

## The pipeline refuses to start, naming `db_root` and `btf_db`

By design. Both database roots must be set, and no script here defaults to a path on any
particular machine. Export `BTF_DB_ROOT` and `BTF_DB`, or pass `--db_root` / `--btf_db`.
See [Installation](Installation).

## `$BTF_DB must end in bioThreatFinder`

Also by design. The build scripts that take `--db-root` append that segment themselves, so
any other name puts the builders and the pipeline in two different trees — and references
under a base the pipeline does not read make every watchlist taxon report `no_reference`,
which reads as "checked and clear" rather than "never built".

## `check_index_coverage.py` exits 3. Is my index broken?

No. **19 of 90 watchlist taxa have no node in the Kraken2 index at all** — the helminths,
several protozoa, most of the dimorphic and opportunistic fungi. That is a property of
PlusPF's library selection, not staleness and not a build error, so rebuilding against the
same libraries changes nothing. Do not spend a day and 123 GB on it.

It is tolerable because all 19 are in the sylph sketch and sylph is on by default: one arm
is blind, not the stage. The check reports separately, and the test suite fails, if a taxon
ever has **no** Kraken2 node *and* no sketched reference — that one would be invisible to
stage 1 entirely.

## A `--check` says my database is stale. What now?

It names what is wrong. Rebuild that database, or regenerate the record if the change was
deliberate. The checks exist because existence is not freshness: a sketch was once 22 days
old, missing sixteen references including every one a divergent filovirus would have
matched, while reporting itself as `ok`.

## `--skip_metabuli false` did nothing

It used to be a silent no-op, and that is worth knowing because it cost this project a
falsification test. On the Nextflow CLI a bare `--flag false` arrives as the **string**
`"false"`, which is truthy in Groovy. `lib/BtfBool.groovy` now coerces these, so the
spelling works — but if you are on an older checkout, that was why.

## The run produced no findings. Is the sample clean?

Not necessarily, and this is the question the tool is most careful about. Check, in order:

1. **The stage-1 coverage table** (`screen/<sample>.screen_summary.tsv`). Did every arm you
   expected actually run? `no_output` means an arm was enabled and produced nothing — and
   for a taxon only that arm would have raised, there would otherwise be no trace.
2. **The unclassified fraction**, in the same table. Sequence whose family is absent from
   the database is not reliably left unclassified: it is confidently assigned to the
   nearest known taxon, at a false-assignment rate approaching 50% on long reads. A novel
   agent appears as familiar background, and this fraction is the only signal that survives
   it.
3. **The database checks** ([Databases](Databases)). A database that does not cover its
   inputs produces exactly this.
4. **`no_reference` rows.** Those are database gaps, not negative results.

## A finding says `supported` and I expected `confirmed_in_silico`

Read the `reasons` column — it names the gate that decided. Three common causes:

- **A marker rule evaluated UNKNOWN.** Not a failure: it means a marker was never
  searched, which is not the same as being absent, so the finding is held one tier down and
  the untested markers are named. 31 of 90 taxa have no rule at all, which caps them by
  construction.
- **An unevaluable vaccine or veto check.** These now cap the tier deliberately: "supported,
  vaccine status untestable" costs a follow-up; "confirmed wild-type poliovirus" costs an
  outbreak response.
- **On ONT, the identity floor.** Every ONT call measured so far failed on
  `identity ... >= 0.9000` and nothing else. See [Validation](Validation) — do not lower it
  on R9.4.1 evidence.

## Can I trust `confirmed_in_silico`?

It means every threshold was met. It does **not** mean a diagnosis: sequence detection does
not establish a viable organism, an active infection, or exposure, and no threshold in the
pack is a measured operating point. It still requires orthogonal laboratory confirmation.

## Why is my taxon `attributed_to_neighbour`?

Because the reads aligned better to a relative in the same competitive index. That is the
tool working: `neighbour_ratio` below ~1 means the relative is the better explanation.
*B. thailandensis* for *B. pseudomallei*, *E. dispar* for *E. histolytica*, and the
B. cereus group for *B. anthracis* are the common ones, and each is a documented trap in the
watchlist's `notes`.

## `-profile docker` fails with `unauthorized`

The images are public, so this normally means an old checkout pointing at images that were
not published yet, or a private registry mirror. Build locally instead:

```bash
docker build -f containers/Dockerfile --build-arg SPEC=sylph -t biothreatfinder-sylph:0.4.0 .
nextflow run . -profile docker --container_prefix biothreatfinder- ...
```

## `--pavian_html` is refused under a container profile

Pavian and two of its dependencies install from GitHub, so no conda spec can declare them
and no image contains them. That process runs in the host environment only. Drop the flag —
the per-sample reports under `pavian/` are still written — or run that step with
`-profile conda`.

## `-profile test` does nothing useful

It exercises the workflow; it does not stand alone, and it cannot be made to by adding a
samplesheet, because stage 1 needs a classifier database and stage 2 the reference pack.
For a test that needs nothing at all: `python3 tests/test_pipeline.py`.

## Bracken did not run on my ONT sample

By design, and the summary says `not_applicable` rather than `not_run`. Bracken's model
assumes one fixed read length, which ONT violates. "Nobody asked for it" and "it cannot
apply here" are different answers to why something was not searched.

## Can I use a GTDB-based Metabuli index?

No, and the pipeline refuses rather than warning. **GTDB taxids are not NCBI taxids**, so
every watchlist lookup would miss — producing not an error but an empty result,
indistinguishable from a clean negative across the entire screen.

## Can I use the GTDB sylph sketch instead of building the pack sketch?

No, for a related reason. sylph reports genome *filenames*, and the taxid is recovered from
the pack's `<taxid>_<name>.fna.gz` naming. GTDB files are named by GTDB accession, so every
lookup would miss and the arm would contribute nothing while the report still said "no
sylph/ANI support".

## Do I have to use the reference lock?

Only if you want the pack the measurements were made on — which you do, if you intend to
compare your results to anything in `docs/`. Without `--from-lock`, each taxon resolves to
whatever NCBI calls latest that day, and two people get two different databases with no
error on either side. Changing one reference has been measured moving read identity
0.9641 → 0.9947.

## Why so many checks?

Because almost every bug found here shares one shape: an empty or degraded result
indistinguishable from a clean negative. A database that opens, an arm that runs, a report
that parses — none of those distinguish a correct answer from a silent subtraction.

## Where do I report a problem?

[`SECURITY.md`](https://github.com/duceppemo/bioThreatFinder/blob/main/SECURITY.md) first,
because the answer depends on what the problem is — and **specimen data must never go into
a public issue**. A wrong result is a correctness bug and belongs in a public issue, where
it is the most useful report this project can receive.
