# Contributing

Issues and pull requests are welcome. Read this first — the conventions here are unusual in
two places, and both are deliberate.

**Before anything else:** do not put specimen data, sample identifiers or a real
investigation's findings into an issue. See [`SECURITY.md`](SECURITY.md).

## The two conventions that are not optional

### 1. A guard must fail loudly, or it is worse than no guard

Most of the bugs found in this project share one shape: *the pipeline produced an empty or
degraded result that was indistinguishable from a clean negative.* An arm that ran and
reported nothing. A reference with no masking result, scored against a mask map that did not
contain it. A sketch that silently omitted sixteen genomes. A marker whose database row was
missing, so the rule evaluated UNKNOWN and the finding settled one tier lower with no
message.

So: when you add a check, make sure a failure is **visible**, and make sure the check
distinguishes *"no"* from *"I could not tell"*. The read-only database checks exit `0` for
yes, `3` for no with the reason named, and `1` for undetermined, and that is not decoration
— conflating the last two is how a missing tool gets reported as a stale database.

The other half of the same rule: **a guard that cries wolf on a healthy system gets
switched off.** Two have been caught doing exactly that and fixed before shipping — one
compared de-duplicated sequence counts and called a correct database stale on 116 of 166
markers, another compared gzip checksums and reported a replaced assembly for byte-identical
sequence fetched by a different route. If your check can fire on a correct installation,
it is not ready.

### 2. Decisions are recorded in `docs/decisions.md`, with the evidence

Every non-obvious choice in this repository has a numbered entry saying what was decided,
what evidence supports it, and what would reverse it. Evidence is tagged
*independent / tool-paper / read-from-source / measured-here / reasoned*, so a reader can
tell a measurement from a judgement. Two entries record decisions that were later reversed
by measurement, and they were left in.

If your change rests on a reason that is not obvious from the diff, add an entry. If it
contradicts an existing entry, say so in the new one and explain what changed. Please do
not quietly delete a finding that turned out to be inconvenient.

## Practicalities

```bash
python3 tests/test_pipeline.py                 # the whole suite, offline, ~30 s
python3 tests/test_pipeline.py --list          # the groups
python3 tests/test_pipeline.py -k markerdb     # one group (substring match, not an expression)
```

- **The suite must pass before and after your change**, and a bug fix should come with the
  check that would have caught it. The suite is fixture-driven, needs no network and no
  databases; cases that compare against a built database skip with a message naming what to
  set, so a bare clone is a valid place to run it.
- **Mutation-test a new check.** Break the thing it is supposed to catch and confirm it
  fails. Three checks in this repo passed a mutation they should have caught, and were
  rewritten; that is cheap to find this way and expensive to find later.
- **Pins live in one place.** A tool version in `modules/local/*.nf` is the source of
  truth; `envs/process/*.yml` and `conf/containers.config` are **generated** from it by
  `scripts/sync_container_specs.py`, and the suite fails if they drift. Change the module,
  then regenerate.
- **No absolute paths.** Nothing may name a path on your machine — not in code, not in a
  comment, not in documentation. `tests/test_pipeline.py -k roots` scans every tracked
  file. Database locations come from `$BTF_DB` and `$BTF_DB_ROOT`, and every script refuses
  rather than defaulting.
- **Nextflow 24.04+**, strict config and DSL syntax. `nextflow run . --help` must compile;
  CI checks it, because a Groovy closure in a local `def` reached a committed file once and
  broke every run.
- **Comments explain why, not what.** The code says what it does. If a line exists because
  of something that was measured, name the measurement or the decision entry.

## Pull requests

One concern per pull request, and a message that says what changed and why rather than
restating the diff. CI runs the suite on two Python versions, compiles `main.nf`, builds
and smoke-tests the container images when a container file changes, and runs ruff and
shellcheck advisorily.

If you are adding a watchlist taxon, a marker, or a reference, note that each of those has a
provenance expectation: a watchlist entry should name its source list, a marker needs a
retrieval route that works and a rule that can evaluate, and a reference should be added to
`packs/*/refs.lock.tsv` so everyone builds the same pack. `docs/roadmap.md` lists the known
gaps if you are looking for somewhere useful to start.

## Conduct

Be straightforward and assume good faith. Disagree about evidence, not about people. That
is the whole code of conduct.
