# Outputs

```
results/
  watchlist/        resolved watchlist + resolution log
  qc/               host_depletion.tsv — how much host got through
  screen/           per-sample stage-1 candidates + screen_summary.tsv
  confirm/          per-(sample,taxon) coverage metrics
  score_evidence/   findings.tsv  <- the machine-readable result
  report/           biothreatfinder_report.html  <- the human-readable one
                    biothreatfinder_summary.tsv
  pavian/           stage-1 reports, converted for Pavian, + optional HTML
  pipeline_info/    timeline, trace, DAG, per-process versions
```

Each process also publishes its own raw output under a lower-case directory named after it
— `kraken2/`, `bracken/`, `sylph_profile/`, `minimap2_competitive/`, `marker_search_aa/`
and so on. **Trimmed FASTQs are not published by default**: they still contain host
sequence, so `--publish_reads` gates them.

## findings.tsv

The result. One row per (sample, taxon), each carrying its tier **and** a `reasons` string
naming the gate that decided it — never a bare verdict. The tiers are described in
[Methods](Methods#evidence-tiers).

Read the gate breakdown, not just the tier. A `supported` call whose only failed gate is
`identity 0.8579 >= 0.9000 (ont/default)` is a very different thing from one that failed
breadth, and the first is what every ONT call in this project has looked like so far.

Columns worth knowing: `n_classifiers` and `classifiers` (how many arms raised it, and
which), `breadth` and `window_occupancy`, `neighbour_ratio` (below ~1 means the reads
prefer a relative), `identity`, `reads_in_masked_regions`, and the marker states per
locus.

## What stage 1 searched, and what it could not

`screen/<sample>.screen_summary.tsv` carries one row per sample per classifier arm, and the
report renders it as **Stage-1 screen coverage**. It answers a question the findings table
cannot: *was this actually searched?*

| status | meaning |
|---|---|
| `ok` | the arm ran and produced a usable report |
| `empty_report` | the arm ran and produced a report with no usable rows |
| `no_output` | the arm was **enabled** and produced nothing |
| `not_applicable` | the arm cannot apply to this sample's platform |
| `not_run` | the arm was disabled for this run |

`not_applicable` is separate from `not_run` on purpose. Bracken is the only arm that takes
it: its model assumes one fixed read length, which ONT violates, so it does not run on long
reads no matter what the operator asked for. Folded into `not_run` it printed "disabled for
this run" on an ONT sample whose kraken2 row read `ok` — inviting the reader to conclude
someone had passed `--skip_kraken2`. "Nobody asked for it" and "it cannot apply here" are
different answers to why something was not searched.

`no_output` is the one that matters. Metabuli, KrakenUniq and sylph each run under
`errorStrategy 'ignore'`, deliberately, so one dead arm cannot sink the run — and the cost
is that a crash leaves a channel padded with the same empty value a *disabled* arm gets.
For a taxon some other arm also raised, the gap shows in `findings.tsv`'s `classifiers`
column. For a taxon **only the dead arm** would have raised — which is the whole reason for
running a divergence-catching arm beside a nucleotide one — there was no trace of it
anywhere. The workflow now tells `parse_screen.py` which arms it enabled, so the difference
is written down and bannered rather than inferred.

Bracken is the exception: it never raises a candidate of its own and exits non-zero by
design on a report with no reads at the requested level — which is what a clean negative
control looks like. Its absence is recorded and not flagged.

`empty_report` is an alarm only for arms that search a comprehensive database. sylph
searches the **pack**, so it answers a sample with no watchlist organism in it with an
empty profile — four of eleven correct answers in one real re-score were bannered as a
classifier that had produced nothing before that was fixed. The run-level version of the
check — an arm empty on *every* sample, which is what a truncated database looks like —
lives in `make_report.py`, the only place that sees every sample at once.

## The unclassified fraction

The same table carries the unclassified fraction per read-based arm and flags it
`HIGH_UNCLASSIFIED` above `--high_unclassified` (0.40). **That number is uncalibrated and
matrix-dependent; set it from your own negatives.**

It is there because of the failure mode
[Perseus](https://doi.org/10.1093/bioinformatics/btag687) measured: sequence whose family
is absent from the database is **not** reliably left unclassified, it is confidently
assigned to the nearest known taxon, at a false-assignment rate approaching 50% on long
reads, sometimes on a median of five k-mers. A novel agent then appears as familiar
background rather than as a gap, and the unclassified fraction is the only signal that
survives it. Triage that fraction — assemble it, protein-search it — rather than reading a
clean findings table as a clean sample.

## Browsing the screen in Pavian

`pavian/` holds every stage-1 arm as a plain 6-column Kraken report, named
`<sample>.<classifier>.pavian.txt`, so [Pavian](https://github.com/fbreitwieser/pavian) can
load a whole run from one directory and line the arms up in its sample-comparison view:

```r
remotes::install_github("fbreitwieser/pavian")
pavian::runApp(port = 5000)   # Data Selection -> browse to <outdir>/pavian
```

These are copies. The export has no output edge into scoring and cannot move a call;
`--skip_pavian_export` turns it off.

It is not merely a convenience, because two arms are not loadable as published:

| arm | as published | why |
|---|---|---|
| kraken2 | **loads** | 8 columns from `--report-minimizer-data`; Pavian's `read_report()` handles it, though `read_report2()` does not |
| bracken | **loads** | already 6-column |
| krakenuniq | **loads** | recognised by its own header |
| metabuli | **rejected** | its report opens with a `#`-prefixed header row, which Pavian reads as a header whose column names it cannot map, so it discards the *entire* file |
| sylph | **unreadable** | not a Kraken report at all: ANI and containment, no tree, no read counts |

The converter drops the k-mer columns, because 6 fields is the only layout both of Pavian's
readers accept. When you want the minimizer evidence, load
`kraken2/<sample>.kraken2.report.txt` directly.

**The sylph export is synthesised, not converted.** Its lineages come from the NCBI taxdump
and its counts are `Sequence_abundance × total_reads`, back-derived so Pavian has something
to draw. They are **not** read assignments and must not be quoted as read counts.
`<sample>.pavian_manifest.tsv` records, per file, whether counts are measured or estimated
and on what basis.

### A static HTML report

`scripts/pavian_report.R` renders Pavian's own standalone template from any run's `pavian/`
directory into one self-contained HTML file that opens with no server and no R:

```bash
# once: Pavian is GitHub-only, so install it where you want it
Rscript -e 'mylib <- "~/R/pavian"; dir.create(mylib, recursive=TRUE)
            .libPaths(c(mylib, .libPaths()))
            remotes::install_github("fbreitwieser/pavian", lib=mylib, upgrade="never")'

scripts/pavian_report.R --indir <outdir>/pavian --lib ~/R/pavian
```

The pipeline can render it too, with `--pavian_html --pavian_r_lib ~/R/pavian`.

**It is off by default and it has to be.** Pavian cannot be installed from conda: it comes
from GitHub, and so do two of its dependencies (`sankeyD3`, `shinyFileTree`). No `conda`
directive can declare it, and installing at run time would need network access and stop
being reproducible. `PAVIAN_REPORT` therefore runs in the **host** environment and is
**refused under any container engine** at launch, rather than failing inside an image with
no R. `PAVIAN_EXPORT` stays unconditional and needs no R, so the reports are always there
whether or not the HTML is.
