# Databases

```bash
scripts/install_databases.sh --route minimal     # ~62 GB
scripts/install_databases.sh --route standard    # + Kraken2 PlusPF, ~177 GB (~700 GB transient)
scripts/install_databases.sh --route full        # + Metabuli, ~286 GB
```

## Downloaded versus built, and why so little is downloaded

| | Downloaded | Built locally |
|---|---|---|
| **what** | the pinned taxonomy snapshot and the curated marker database — **166 MB of archive**, [doi:10.5281/zenodo.23602998](https://doi.org/10.5281/zenodo.23602998) | the 43 GB reference pack, the host genome, the Kraken2 / Metabuli indexes |
| **why that split** | NCBI overwrites `new_taxdump.tar.gz` **in place**, so the pinned snapshot is gone upstream — and every taxid the pipeline reports keys to it. The marker database is curation, not compute | `packs/<pack>/refs.lock.tsv` pins every accession **and** the checksum of its decompressed sequence, so a clone builds the same pack and can prove it. A lock is auditable; a 115 GB archive can only be taken on trust |
| **verified by** | sha256, before anything is unpacked | the lock's per-reference checksum, then the read-only checks below |

Each archive carries `DISCLAIMER.md` and a `BTF_PROVENANCE.txt` inside, because an archive
outlives the page it was downloaded from. `--from-dir DIR` installs the same archives from
a local copy, verifying the same checksums — so a copy from a colleague is as trustworthy
as a download.

**Kraken2 is built, not downloaded, and the upstream prebuilt indexes are not
equivalent.** The capped variants (`k2_pluspf_16_GB` and friends) are built with
`--max-db-size` and lose sensitivity, and any upstream index carries its own taxonomy
snapshot rather than the pinned one everything else here keys to.

## The reference lock

```bash
python3 scripts/build_pack_refs.py --pack human --from-lock packs/human/refs.lock.tsv
```

Without `--from-lock`, each watchlist taxon resolves to **whatever NCBI calls latest that
day**. The selection rule is deterministic; the catalogue is not. So two people build two
different reference packs from the same watchlist, **with no error on either side** — and
then the sylph sketch covers different genomes, stage 2 aligns against different
assemblies, and every identity, breadth and limit-of-detection figure in `docs/` describes
a deployment nobody else has.

The size of that risk was measured: changing **one** reference moved read identity
0.9641 → 0.9947 and best-reference breadth 0.877 → 0.982, which is the difference between
a call and a miss.

The lock holds 375 assemblies across 237 taxa, each pinned by accession and by the md5 of
its **decompressed** sequence. The decompressed checksum matters: the gzip container
differs between download routes for byte-identical sequence, so verifying the container
would report a replaced assembly to anyone who fetched it a different way — and a guard
that cries wolf on a healthy database gets switched off.

Four behaviours, each verified against NCBI rather than reasoned:

- the locked accession is fetched, **not** the top-ranked one, and exactly as many as the
  lock names — `--max-per-taxon` must not trim a locked taxon;
- a checksum mismatch fails the build, naming the accession;
- an accession NCBI no longer serves fails the build with `THIS PACK IS NOT THE LOCKED
  PACK`, after doing everything it can;
- a taxon locked as having **no assembly anywhere** is recorded as such — four
  near-neighbours in the human pack exist nowhere, and that is a different fact from "the
  lock predates this taxon".

Regenerate after deliberately changing the pack:

```bash
python3 scripts/build_pack_refs.py --pack human --write-lock packs/human/refs.lock.tsv --dry-run
```

## What each artifact is

| Artifact | Size | Needed for |
|---|---|---|
| `taxdump/2026-09-29/` | 2.2 GB | watchlist resolution, every classifier's taxonomy. 3,016,439 nodes. **Published**, because it cannot be re-fetched |
| `packs/<pack>/refs/` | 43 GB | stage 2. 375 references with sequence, 90 competitive bundles + 180 minimap2 indexes |
| `packs/<pack>/refs/masked/` | — | the breadth denominators. 20.0 Mbp masked across all 375 |
| `packs/<pack>/markers/` | 61 MB | the marker tier. 166 rows, 146 built, 33 profile HMMs. **Published** |
| `sylph/<pack>_pack_refs.syldb` | 177 MB | the sylph arm. 375 genomes, 6 s to rebuild |
| `host/<pack>/` | 17 GB | host depletion, and masking's host pass |
| `kraken2/pluspf/` | 115 GB | the primary stage-1 arm. `hash.k2d` is all but 0.2 GB of it; **~700 GB free is needed to build it**, and the finished index is what remains |
| `metabuli/refseq` | 110 GB | the opt-in DNA+AA arm |

The sylph arm uses the **pack's own** sketch, not a GTDB sketch. They are not
interchangeable: sylph reports genome *filenames*, and `parse_screen.py` recovers a taxon
by reading the taxid out of the pack's `<taxid>_<name>.fna.gz` naming. GTDB files are named
by GTDB accession, so every lookup would miss and the arm would contribute nothing while
the report still said "no sylph/ANI support".

## Is a database still current?

Existence is not freshness. A database is built once and then goes on existing while the
pack it describes keeps growing — which is how a sketch came to be 22 days stale, missing
sixteen references including every one a divergent filovirus would have matched, while
reporting itself as `ok`.

Each of these writes nothing. **0** means yes, **3** means no and it names what is wrong,
**1** means it could not tell — and conflating the last two is how a missing tool gets
reported as a stale database.

```bash
python3 scripts/build_marker_db.py      --pack human --check    # covers markers.tsv, on every claimed tier
bash    scripts/build_sylph_db.sh       --pack human --check    # covers every reference the MANIFEST names
python3 scripts/check_pack_bundles.py   --pack human            # do the competitive bundles still match the manifest
python3 scripts/check_index_coverage.py --pack human            # can Kraken2 report every watchlist taxon at all
python3 scripts/check_metabuli_coverage.py --pack human         # the same question for the Metabuli arm
python3 scripts/mask_pack_refs.py       --pack human --dry-run  # a result for every reference, and a reference for every result
python3 scripts/check_threshold_coherence.py                    # which per-taxon gates can never bind
```

Three of them also run **at launch** — the only consumer of these databases is a run, so a
run is where the operator is told. They warn and never refuse, because a stale database
still answers and the operator may know why.

**`check_index_coverage.py` exits 3 on the human pack by design.** 19 of 90 watchlist taxa
have no node in the Kraken2 index at all: the helminths, several protozoa, and most of the
dimorphic and opportunistic fungi. That is **not** staleness and not a build error — it
follows from PlusPF's library selection, so rebuilding against the same libraries changes
nothing. It is tolerable only because all 19 are in the sylph sketch and sylph is on by
default: **one arm is blind, not the stage.** A taxon with no Kraken2 node *and* no sketched
reference would be invisible to stage 1 entirely; the check reports that case separately and
the test suite fails if it ever appears.

**The masking check asks its question both ways round (D62).** "Does every reference have a
result" was the only direction for a long time, and the other one caught 31 mask BEDs on
2026-10-26 that belonged to references the manifest no longer named — left behind because
the manifest merge carries every earlier row forward, which is the right behaviour for an
interrupted run and the wrong one for a reverted reference set. They made the artifact
disagree with its own documentation: 20,148,007 bp over 406 results, against the 20,040,354
over 375 that `docs/` correctly reports. Whether an orphan is inert is measured rather than
assumed, because `coverage_stats.py --mask-bed` globs the whole directory into one
contig→intervals map: an orphan whose contigs are still in the pack masks live sequence no
manifest row accounts for, and that exits non-zero, while a dead one is a warning.
`--prune-orphan-masks` removes them.

Only the sylph check needs a tool on `PATH` (`conda activate btf-sylph`); the rest are
stdlib Python and run from a bare clone. That is deliberate — a check that needs a
configured deployment is a check nobody runs on a fresh clone.

**Metabuli's answer is 3 as well**, and more sharply: 63 of 90 watchlist taxa have no
sequence at or below them in the prebuilt database that ships — every fungus, 19 of 20
parasites. The arm is off by default, so this is not a run-time defect; it is what the arm
would answer if switched on for those taxa. A different Metabuli database would answer
differently.

**The bundle check reads the bundles, and that is the point.** Its staleness predicate is
*imported* from `build_pack_refs.py`, so there is one predicate rather than two that can
drift. But only 10 of the 90 bundles carry the `.inputs.tsv` record that predicate prefers —
the other 80 were built before it existed, and for those it can only compare timestamps,
which a restore from backup or a re-download of the same accession defeats. So the default
opens every bundle and compares each taxid prefix's contig count and base total against the
manifest: no sidecar needed, and stronger than one, because the sidecar records what the
builder *believed* it wrote. 73 s for the pack, 2.1 GB read; `--quick` gives the
metadata-only answer, which on this pack is honestly **1**.

That pass also re-applies the two checks that otherwise only ever run during a rebuild —
identical sequence under two taxid prefixes, which would give every read over it MAPQ 0, and
a header with no `taxid|` prefix, which no read can be attributed to. All 90 bundles pass
both. What it does **not** answer is whether the manifest matches `refs.lock.tsv`; the lock
outranks the disk, and that comparison needs `build_pack_refs.py --from-lock`, which reaches
NCBI.

## Publishing a database set

```bash
scripts/package_databases.sh --out-dir DIR
```

Builds the two publishable archives, refuses to package a marker database that fails its
own freshness check, computes the checksums into `assets/databases.manifest.tsv`, and
points at `assets/zenodo_record.md` for the deposit metadata. The operational detail —
sizes, rebuild costs, what is re-creatable and what is not — is in
[`docs/operations.md`](https://github.com/duceppemo/bioThreatFinder/blob/main/docs/operations.md).
