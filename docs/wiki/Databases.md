# Databases

```bash
scripts/install_databases.sh --route minimal     # ~62 GB
scripts/install_databases.sh --route standard    # + Kraken2 PlusPF, ~323 GB
scripts/install_databases.sh --route full        # + Metabuli, ~433 GB
```

## Downloaded versus built, and why so little is downloaded

| | Downloaded | Built locally |
|---|---|---|
| **what** | the pinned taxonomy snapshot and the curated marker database — **166 MB of archive**, [doi:10.5281/zenodo.23602998](https://doi.org/10.5281/zenodo.23602998) | the 43 GB reference pack, the host genome, the Kraken2 / Metabuli indexes |
| **why that split** | NCBI overwrites `new_taxdump.tar.gz` **in place**, so the pinned snapshot is gone upstream — and every taxid the pipeline reports keys to it. The marker database is curation, not compute | `packs/<pack>/refs.lock.tsv` pins every accession **and** the checksum of its decompressed sequence, so a clone builds the same pack and can prove it. A lock is auditable; a 261 GB archive can only be taken on trust |
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
| `kraken2/pluspf/` | 144 GB | the primary stage-1 arm |
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
python3 scripts/build_marker_db.py    --pack human --check      # covers markers.tsv, on every claimed tier
bash    scripts/build_sylph_db.sh     --pack human --check      # covers every pack reference
python3 scripts/check_index_coverage.py --pack human            # can the index report every watchlist taxon at all
python3 scripts/mask_pack_refs.py     --pack human --dry-run    # does every reference have a masking result
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

Only the sylph check needs a tool on `PATH` (`conda activate btf-sylph`); the rest are
stdlib Python and run from a bare clone. That is deliberate — a check that needs a
configured deployment is a check nobody runs on a fresh clone.

**Two gaps, both known.** Metabuli has no freshness or coverage check at all; it is opt-in
and off by default. And the competitive bundles' predicate is a rebuild decision inside
`build_pack_refs.py` rather than a standalone check, so asking "are my bundles current?"
means invoking the builder.

## Publishing a database set

```bash
scripts/package_databases.sh --out-dir DIR
```

Builds the two publishable archives, refuses to package a marker database that fails its
own freshness check, computes the checksums into `assets/databases.manifest.tsv`, and
points at `assets/zenodo_record.md` for the deposit metadata. The operational detail —
sizes, rebuild costs, what is re-creatable and what is not — is in
[`docs/operations.md`](https://github.com/duceppemo/bioThreatFinder/blob/main/docs/operations.md).
