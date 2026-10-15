#!/usr/bin/env Rscript
# Render Pavian's own static HTML report from a run's pavian/ directory.
#
# WHY THIS IS A SCRIPT AND NOT A PIPELINE PROCESS
# -----------------------------------------------
# Pavian is not installable from conda. It comes from GitHub, and two of its
# dependencies (sankeyD3, shinyFileTree) are GitHub-only as well, so a Nextflow
# `conda` directive cannot declare it and a process that installed it at run
# time would need network access and would not be reproducible. The pipeline
# therefore emits the REPORTS (always, no R required) and leaves the HTML to
# this script, which you point at any run afterwards.
#
# The HTML is Pavian's, not ours: it renders inst/pavian-report.Rmd from the
# installed package, which is the same template the Shiny app's "Generate
# report" button uses. Nothing here reimplements a view.
#
# Usage:
#   scripts/pavian_report.R --indir <outdir>/pavian [--out report.html]
#                           [--title ...] [--author ...] [--lib <R library>]
#                           [--no-sankey]
#
# First time, install Pavian into a library of your choice:
#   Rscript -e 'mylib <- "~/R/pavian"; dir.create(mylib, recursive=TRUE)
#               .libPaths(c(mylib, .libPaths()))
#               remotes::install_github("fbreitwieser/pavian", lib=mylib,
#                                       upgrade="never")'
# then pass --lib ~/R/pavian (or leave it out if it is on the default path).

args <- commandArgs(trailingOnly = TRUE)

getopt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(flag, " needs a value")
  args[i + 1]
}
has_flag <- function(flag) flag %in% args

indir  <- getopt("--indir")
if (is.null(indir)) stop("--indir <pavian directory> is required")
indir  <- normalizePath(indir, mustWork = TRUE)
outfile <- getopt("--out", file.path(indir, "pavian_report.html"))
title   <- getopt("--title", paste("bioThreatFinder - stage-1 classification -",
                                   basename(dirname(indir))))
author  <- getopt("--author", "bioThreatFinder")
extra_lib <- getopt("--lib")

if (!is.null(extra_lib)) .libPaths(c(normalizePath(extra_lib, mustWork = TRUE), .libPaths()))

if (!requireNamespace("pavian", quietly = TRUE)) {
  stop("pavian is not installed or not on .libPaths(). See the header of this ",
       "script for the one-line install, and pass --lib if you installed it ",
       "into a private library.")
}

# Only the converted reports. The manifests and README live in the same
# directory and Pavian rejects them as malformed anyway, but naming the files
# explicitly keeps that out of the log and lets us set readable sample names.
report_files <- sort(Sys.glob(file.path(indir, "*.pavian.txt")))
if (length(report_files) == 0)
  stop("no *.pavian.txt files in ", indir,
       " - run the pipeline (PAVIAN_EXPORT) or bin/to_pavian.py first")

# "case_1-1.kraken2.pavian.txt" -> "case_1-1 (kraken2)", so the report's tables
# are readable instead of showing full filenames.
base  <- sub("\\.pavian\\.txt$", "", basename(report_files))
samp  <- sub("\\.[^.]+$", "", base)
clf   <- sub("^.*\\.", "", base)
names_nice <- paste0(samp, " (", clf, ")")

# sample_data.csv is read by Pavian's own directory loader (with sep=";", despite
# the name), so writing it here also gives these readable names to the Shiny app
# when someone browses to this directory later.
sample_data <- data.frame(
  FormatOK = TRUE, Include = TRUE,
  Name = names_nice,
  ReportFile = basename(report_files),
  ReportFilePath = report_files,
  stringsAsFactors = FALSE
)
utils::write.table(sample_data[, c("Name", "ReportFile")],
                   file.path(indir, "sample_data.csv"),
                   sep = ";", row.names = FALSE, quote = FALSE)

cat("[pavian] loading", length(report_files), "report(s)\n")
reports <- pavian::read_reports(report_files, names_nice)

rmd <- system.file("pavian-report.Rmd", package = "pavian")
if (!nzchar(rmd)) stop("pavian is installed but pavian-report.Rmd is missing from it")

# Render in a temp copy, exactly as the Shiny app does, so the package
# directory is never written to.
tmp_rmd <- file.path(tempdir(), "pavian-report.Rmd")
file.copy(rmd, tmp_rmd, overwrite = TRUE)

outfile <- ifelse(grepl("^/", outfile), outfile, file.path(getwd(), outfile))
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)

rmarkdown::render(
  tmp_rmd,
  output_file = outfile,
  params = list(
    doc_title = title,
    doc_author = author,
    doc_date = date(),
    set_name = basename(dirname(indir)),
    all_data_loaded = TRUE,
    sample_data = sample_data,
    reports = reports,
    include_sankey = !has_flag("--no-sankey"),
    filter_taxa = NA
  ),
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)

cat("[pavian] wrote", outfile, "\n")
