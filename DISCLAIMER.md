# Disclaimer

**bioThreatFinder is research software under active development. It is not a validated
diagnostic assay, and no result it produces is a diagnosis, a clearance, or a basis for
clinical, regulatory, trade or biosecurity action on its own.**

This file ships with the software and with every published database archive, because the
databases carry the same limits as the code that reads them.

## What is not established

- **No threshold is a measured operating point.** Every number in `packs/*/pack.yaml` is a
  reasoned default drawn from the literature, not a value calibrated against characterised
  material. Three are explicitly marked PROVISIONAL in the code, and
  `scripts/check_threshold_coherence.py` reports 75 of 89 taxa carrying a `min_unique_reads`
  gate that can never bind.
- **Sensitivity and specificity are not characterised for any intended use.** What has been
  measured is recorded in `docs/verification.md` and `docs/benchmarking.md`, including a
  limit of detection on one organism in one simulated matrix and a run against the
  ZymoBIOMICS standards. That is a long way from a validated performance claim.
- **The marker-gene tier addresses gaps for which no published method exists** — vaccine
  versus wild-type calling from metagenomic data among them. That is why it was built, and
  it is also why it is unvalidated. It must not be presented as established practice.
- **The watchlist is not authoritative.** It derives from published notifiable and
  high-consequence lists but has not been reconciled against the legal source for any
  jurisdiction. `docs/roadmap.md` records what reconciliation requires. Risk-group
  assignments in particular must come from the relevant national registry, not from here.
- **Coverage is incomplete and partly measured.** 19 of 90 watchlist taxa have no node in
  the Kraken2 index and are reachable only through the sylph arm (D48). 31 of 90 have no
  marker rule, so the top evidence tier is unreachable for them by construction. Four near
  neighbours have no reference genome anywhere, which weakens discrimination for their
  targets. Each of these is recorded rather than smoothed over, and each affects what an
  absent result means.
- **A negative is weaker than it looks.** The failure mode this project keeps finding in
  itself is a result that is empty or degraded and *indistinguishable from a clean
  negative*. Several guards exist specifically to make that visible; none of them can make
  absence of evidence into evidence of absence.

## What the databases are and are not

The published archives are **derived data**: a dated snapshot of NCBI taxonomy, and
reference and marker sequences retrieved from public NCBI records. They are redistributed
to make a specific analysis reproducible — not as an authoritative or curated reference
resource, and not as a substitute for the primary records, which remain the citable source.
Sequence records are reproduced under NCBI's terms; where an individual submitter asserts
additional conditions, those apply to that record.

A snapshot is a snapshot: taxonomy, assemblies and gene annotations all change, and the
pinned taxonomy is deliberately frozen so that results remain comparable. It is therefore
out of date by design, and progressively more so.

## Appropriate use

Intended for method development, benchmarking and research, by users competent to judge the
output — which includes reading the per-gate evidence rather than the tier label. Any
result bearing on a real decision requires confirmation by a validated, fit-for-purpose
method, and a human assessment of the evidence behind the call.

Nothing here should be used to decide that a sample is safe.

## No warranty

Provided as is, without warranty of any kind, express or implied, under the terms in
[`LICENSE`](LICENSE). The authors and the Canadian Food Inspection Agency accept no
liability for any decision taken on the basis of its output. Views and outputs are those of
the software and its authors, and are not an official position of the Canadian Food
Inspection Agency or the Government of Canada.
