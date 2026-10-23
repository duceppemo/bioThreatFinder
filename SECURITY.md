# Reporting problems

## Do not put specimen data in a public issue

This matters more here than the usual software-security advice, so it comes first.

**Never attach or paste sequence reads, sample identifiers, patient or premises
identifiers, case metadata, or the findings table from a real investigation into a GitHub
issue, pull request or discussion.** A public tracker is a public record, metagenomic reads
from a clinical or veterinary specimen can carry host sequence that identifies the host,
and a `findings.tsv` naming a high-consequence organism alongside a date and a site is
information about an incident, not a bug report.

To report a bug that only shows up on real data:

1. Reproduce it with simulated or public data — `scripts/simulate_benchmark.py`, or a
   ZymoBIOMICS accession — and report that.
2. If it cannot be reproduced that way, open an issue describing **only** the shape of the
   problem: the process, the message, the pipeline version, and what the inputs were like
   (platform, read count, rough composition). Say that you have data you cannot share.
3. If the problem is in output that cannot be described without the data, email
   <marc-olivier.duceppe@inspection.gc.ca> and we will arrange a channel. Expect to be
   asked whether the data can be shared at all.

Redact before you send anything, and assume redaction is imperfect: sample names leak into
file names, read headers, and the HTML report.

## Software vulnerabilities

This is a research pipeline, not a network service. It has no server, no authentication, no
user accounts, and it opens no ports. The realistic attack surface is:

- **Input parsing.** FASTQ, samplesheet CSV, `pack.yaml`, `markers.tsv`, `marker_rules.yaml`
  and the classifier reports are all parsed by code in `bin/` and `scripts/`. A crafted
  input that causes command injection, a path traversal, or arbitrary code execution is a
  genuine vulnerability — report it.
- **Database download and unpack.** `scripts/install_databases.sh` fetches archives and
  unpacks them. Checksums are verified before unpacking; a way to get past that, or a path
  traversal out of the unpack directory, is a vulnerability.
- **Dependencies.** Report a known-vulnerable pinned dependency as an ordinary issue; the
  pins are in `modules/local/*.nf` and `envs/`.

Report any of the above **privately first**, by email to
<marc-olivier.duceppe@inspection.gc.ca>, rather than in a public issue. Include what you
ran, what happened, and the version. There is no bounty and no formal SLA; expect an
acknowledgement within a couple of weeks, and say in the subject line if you believe the
issue is time-sensitive.

## What is not a vulnerability

- **The pipeline reaching a wrong conclusion.** That is a correctness bug, and it belongs
  in a public issue where it can be discussed and measured. It is also the most useful kind
  of report this project can receive — see `docs/verification.md` for what is and is not
  checked, and `DISCLAIMER.md` for why no output here is a diagnosis.
- **A threshold being wrong.** No threshold in `packs/*/pack.yaml` is a measured operating
  point; they are reasoned defaults, documented as such. Measurements that would improve
  them are very welcome.
- **The watchlist being incomplete or out of date.** Known and recorded. It is not an
  authoritative legal source for any jurisdiction.

## Supported versions

The `main` branch is what is supported. There are no backports.
