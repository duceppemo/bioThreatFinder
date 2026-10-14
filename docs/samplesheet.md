# The samplesheet

The one file every run needs, and until 2026-10-13 the only part of the interface with
no description outside `main.nf`'s error strings. `docs/validation.md` has recorded the
gap for longer than that.

```csv
sample,platform,fastq_1,fastq_2,role
SPEC001,illumina,/data/SPEC001_R1.fastq.gz,/data/SPEC001_R2.fastq.gz,sample
NTC01,illumina,/data/NTC01_R1.fastq.gz,/data/NTC01_R2.fastq.gz,control
ONT042,ont,/data/ONT042.fastq.gz,,sample
```

Passed as `--input`. Parsed by `parseSamplesheet` in `main.nf`; also read directly by
`score_evidence.py` and `make_report.py`, which is why it is authoritative rather than
advisory (see `role`).

## Columns

| column | required | values | notes |
|---|---|---|---|
| `sample` | yes | any non-empty string | Unique. Becomes `meta.id` **verbatim** and is used as a filename component (`<sample>.<taxid>.coverage.tsv`), so avoid `/`, whitespace and shell metacharacters. It is carried unaltered into every output; nothing renames, suffixes or sanitises it. |
| `platform` | yes | `illumina` \| `ont` | Case-insensitive. Selects the read-identity floor, the alignment-length fraction and the simulator-independent gates. A wrong value is not an error and will silently judge the sample against the other platform's thresholds. |
| `fastq_1` | yes | path | Must exist at launch. |
| `fastq_2` | no | path or empty | Empty means single-end. **Must** be empty for `ont`, which is an error rather than a warning. |
| `role` | no | `sample` \| `control` | Defaults to `sample`. Anything else is a hard error. |

Unknown columns are **ignored silently** — `parseSamplesheet` checks only that
`sample`, `platform` and `fastq_1` are present. A misspelled header (`Role`, `platfrom`)
therefore does not fail the run; it falls back to the default. Check your header.

## `role`, and why it is authoritative

`role=control` marks a negative control. The samplesheet is the **only** statement of
this that the pipeline trusts, because counting controls from pipeline output answers a
different question and answers it wrongly in the best case: a spotless NTC raises
nothing, produces no rows anywhere, and would count as zero controls — putting
"NO NEGATIVE CONTROL IN RUN" on every finding of the cleanest possible run.

Two things follow from a row being a control:

1. **Its content is never an actionable finding.** Any tier that asserts presence is
   reclassified to `control_background`. A control's content is QC about the run, not a
   finding about a specimen.
2. **It suppresses.** For each taxon, the highest control RPM in the run sets a floor:
   a sample is tiered `suppressed_control` when its RPM does not exceed
   `control_ratio` (default 5.0) times that value.

### The scope of a control is the whole run

There is no way to say "this control covers only these samples". The highest control
value for a taxon applies to **every** sample in the run, across every scenario, plate
or extraction batch in that samplesheet.

This is deliberate and it is load-bearing — see `docs/decisions.md` D21. It is also the
sharpest edge in the file, so know these consequences:

- A control that genuinely contains an organism suppresses that organism **everywhere**
  in the run. If one extraction batch's blank is contaminated, every other batch in the
  same samplesheet inherits its floor.
- The failure arrives when the operator does the right thing: adding a control
  simultaneously removes the "no negative control" caveat and can suppress a true
  positive.
- A 5x margin means a control at 1,000 RPM owns everything below 5,000 RPM — a sample
  can carry **more** of the organism than the control does and still be suppressed
  (D22).

**None of this is silent any more.** Every run prints a `CONTROL SUBTRACTION AUDIT`
table to `score.log` naming each taxon any control raised, the control responsible, and
the implied floor — including taxa that suppressed nothing yet, marked `LATENT`. Each
finding carries a `control_subtraction` column (`cleared:` / `suppressed:` /
`exceeds_control:` / `not_evaluated:`), and the HTML report has a matching section.
**Read that table before trusting a negative.**

If one control must not speak for the whole run, split the samplesheet and run the
groups separately. A per-group scoping column was designed and rejected — the
measurement is in D21 and the replacement belongs with roadmap item B.

## Benchmark samplesheets

`scripts/simulate_benchmark.py` writes a runnable samplesheet beside its manifest. It
refuses to build a batch in which a scenario declaring `index_hop_fraction` — a
deliberately contaminated control — shares the batch with any other scenario, because
that control would impose a floor on its target across all of them. Run it alone with
`--scenario`.
