# Zenodo deposit metadata

Ready to paste. `scripts/package_databases.sh` builds the archives this describes and
prints a pointer here; this file is the single source of truth for the record's wording, so
edit it rather than the script.

**One record, both files.** They are installed together, they version together, and 166 MB
does not need splitting. When a second domain pack arrives, add its marker archive as a
*new version of this same record* rather than a new record — Zenodo keeps a concept DOI
that always resolves to the latest version, plus a DOI per version. Cite the concept DOI in
the README and `CITATION.cff`; put the **version** DOI's file URLs in
`assets/databases.manifest.tsv`, because the installer needs exact files with the checksums
it records, not "whatever is newest".

---

## Resource type

**Dataset.** Not "Software" — the software is the GitHub repository, and linking the two is
what `Related identifiers` below is for.

## Title

> bioThreatFinder databases: a pinned NCBI taxonomy snapshot (2026-09-29) and a curated
> marker-gene database for high-consequence pathogen screening

Why this shape: the snapshot date is in the title because it is the single fact that makes
the deposit necessary, and someone looking for a dated NCBI taxonomy snapshot should be
able to find it without reading the description.

## Authors

| Field | Value |
|---|---|
| Name | Duceppe, Marc-Olivier |
| Affiliation | Canadian Food Inspection Agency |
| ORCID | `0000-0003-2130-0427` — in `CITATION.cff` and `.zenodo.json`, so every future release carries it |

## Version

`2026-10-23`

Date-based on purpose. The content is a snapshot, so a semantic version would imply an API
that these files do not have. The software version belongs in the description, not here.

## Publication date

The deposit date.

## Language

English (eng)

## License

**Creative Commons Attribution 4.0 International (CC-BY-4.0).**

The repository's MIT licence covers the code, not this. CC-BY is the conventional choice
for a derived-data deposit and keeps attribution attached to the curation, which is the part
that is actually yours: the marker catalogue, the retrieval decisions per marker, and the
record of which markers cannot be built and why. CC0 is the alternative if you would rather
impose nothing at all; the underlying NCBI records carry no such requirement either way.

## Keywords

```
metagenomics
biosurveillance
pathogen detection
high-consequence pathogens
reference database
NCBI taxonomy
taxonomy snapshot
marker genes
vaccine strain discrimination
Nextflow
bioinformatics pipeline
reproducibility
public health microbiology
one health
biothreat
```

The first half is what a domain reader searches for; the second half is what someone
looking for *a reproducible pipeline's pinned databases* searches for, which is a different
person with a different query.

## Dates

Add a date of type **Collected**: `2026-09-29` — the day the taxonomy snapshot was taken.
Zenodo's `Publication date` will be the deposit date, which is weeks later and says nothing
about what is in the files.

## Related identifiers

| Relation | Identifier | Type |
|---|---|---|
| Is supplement to | `https://github.com/duceppemo/bioThreatFinder` | Software |
| Is derived from | `https://www.ncbi.nlm.nih.gov/taxonomy` | Dataset |
| Is derived from | `https://www.ncbi.nlm.nih.gov/datasets/` | Dataset |
| Is required by | *(the software's own Zenodo DOI, once the GitHub→Zenodo integration is switched on and a release is tagged)* | Software |

That last one is worth doing in the same sitting: enabling the integration and tagging a
release gives the *code* a DOI too, and then the pair cite each other. A dataset whose
software is uncitable is half a record.

## Description

Paste as-is. It repeats the disclaimer deliberately: a record page is read by people who
will never open the archive, and the archive is opened by people who will never see the
record page.

---

bioThreatFinder is a Nextflow pipeline for metagenomic screening and in-silico confirmation
of high-consequence pathogens across human, animal and plant health. This deposit holds the
two database artifacts that **cannot be re-created from public sources**. Everything else
the pipeline needs is rebuilt by its own scripts.

**What is here**

1. `btf-taxdump-2026-09-29.tar.gz` — the NCBI `new_taxdump` snapshot of 2026-09-29,
   3,016,439 nodes. NCBI overwrites `new_taxdump.tar.gz` in place, so this exact snapshot
   can no longer be fetched from NCBI. Every taxid the pipeline reports keys to it: the
   watchlist is resolved against it, the Kraken2 index is built against it, and a different
   snapshot silently moves clade membership — a copy 6.5 months older was measured at
   +288,764 taxids.
2. `btf-markers-human-2026-10-21.tar.gz` — the curated marker-gene database for the human
   domain pack: 166 markers across nucleotide, protein and profile-HMM evidence tiers, with
   33 profile HMMs. `markers.manifest.tsv` inside is authoritative for which markers were
   actually built; a marker absent from it evaluates UNKNOWN and never promotes a finding.
   16 of the 166 can only ever be UNKNOWN, and the reason for each is recorded rather than
   hidden.

**What is deliberately not here.** The 43 GB reference pack is not deposited. It is pinned
instead by `packs/human/refs.lock.tsv` in the repository, which names every assembly
accession and the md5 of its *decompressed* sequence, so a clone rebuilds the exact pack
these measurements were made on and can prove it. A lock is auditable; a multi-gigabyte
archive can only be taken on trust. The Kraken2 (115 GB) and Metabuli (110 GB) indexes are
built from their own scripts for the same reason.

**Use**

```
export BTF_DB=/fast/disk/bioThreatFinder
scripts/install_databases.sh --route minimal
```

The installer verifies each archive's sha256 **before** unpacking it, because the failure it
prevents is quiet: a truncated taxonomy does not crash the watchlist resolver, it resolves
fewer taxa and reports success. The checksums are recorded in
`assets/databases.manifest.tsv` in the repository.

**Provenance.** Each archive contains a `BTF_PROVENANCE.txt` naming the snapshot date and
the pipeline commit that packaged it, and a copy of `DISCLAIMER.md`.

**Disclaimer.** bioThreatFinder is research software under active development, **not a
validated diagnostic assay**. No threshold in it is a measured operating point; sensitivity
and specificity are not characterised for any intended use; the watchlist derives from
published notifiable and high-consequence lists but is **not** an authoritative legal
source for any jurisdiction; and the marker tier addresses gaps for which no published
method exists, which is precisely why it is unvalidated. Coverage is incomplete and partly
measured: 19 of 90 watchlist taxa are reachable through one stage-1 arm only, and 31 of 90
cannot reach the top evidence tier at all. These archives are **derived data** — a dated
snapshot of public NCBI records, redistributed for reproducibility rather than as a curated
reference resource — and are out of date by design. **No output of this software should be
used to decide that a sample is safe.** Any result bearing on a real decision requires
confirmation by a validated, fit-for-purpose method. See `DISCLAIMER.md`, included in both
archives.

Software, documentation and the decision record:
https://github.com/duceppemo/bioThreatFinder

---

## Additional notes (optional field)

> Produced for bioThreatFinder <TAG>, commit <SHA>. File checksums are recorded in
> `assets/databases.manifest.tsv` in the repository; Zenodo's own md5 is shown per file on
> this page. The taxonomy snapshot is pack-independent; the marker database is specific to
> the human domain pack.

## Funding

Add the CFIA project or grant if this work sits under one — it is the field most often left
empty and most often wanted later by whoever has to report on outputs.

---

## After depositing

1. Copy each file's download URL — Zenodo serves them as
   `https://zenodo.org/records/<record-id>/files/<filename>?download=1` — into the `url`
   column of `assets/databases.manifest.tsv`. Use the **version** record's id, not the
   concept DOI: the installer verifies the sha256 it has recorded, and "latest" would
   eventually not match.
2. Add the **concept** DOI to `CITATION.cff` and to the README's install section.
3. Commit both. `scripts/install_databases.sh` refuses to download while the manifest says
   `PENDING`, so this is the step that switches the public route on.
