# Session handoff — 2026-09-29

Written because the session was ending with long builds in flight. Read this first
on resume.

## Jobs running detached (they survive the session dying)

Verified `PPID 1`, so they are not children of the shell. Check progress with:

```bash
ls -t /media/8tb/bioThreatFinder/logs/ | head
tail -f /media/8tb/bioThreatFinder/logs/build_kraken2_<latest>.log
```

| Job | State when handed off | How to tell it finished |
|---|---|---|
| **Kraken2 PlusPF build** | downloading libraries (~640 MB in; bacteria is the long pole) | log ends with the Bracken section and a SUMMARY; `hash.k2d` exists |
| **Marker DB full pass** | recovering the manifest (see below) | `markers.manifest.tsv` has **163** rows |

If either died, both are resumable and skip completed work:

```bash
export PATH=/home/bioinfo/miniconda3/envs/btf-kraken2/bin:$PATH
bash scripts/build_kraken2_db.sh --db-root /media/8tb \
     --taxdump /media/8tb/bioThreatFinder/taxdump/2026-09-29 --threads 40

export PATH=/home/bioinfo/miniconda3/envs/btf-refs/bin:$PATH
python3 scripts/build_marker_db.py --pack human --db-root /media/8tb --threads 24
```

## Environments — and a trap to avoid

| Env | Holds | Note |
|---|---|---|
| `btf-refs` | datasets, minimap2, seqkit, samtools, mmseqs2, hmmer, mafft | all current **except** its kraken2, which is 2.0.8-beta — do not use it for kraken2 |
| `btf-kraken2` | **kraken2 2.17.2**, bracken | deliberately separate |

**Why separate:** installing kraken2 into `btf-refs` alongside mmseqs2/hmmer/mafft made
the solver resolve **kraken2 2.0.8-beta (2019)**, which has no
`--report-minimizer-data` — the flag the whole stage-1 design rests on. An index built
with it would work and carry none of the evidence the pipeline scores on. Audited the
rest of `btf-refs` afterwards: everything else was current, so this was specific to
kraken2's perl/wget dependencies. `build_kraken2_db.sh` now refuses anything below 2.1.0
and verifies the flag is actually advertised.

Also: `btf-refs`' wget is broken (linked against `libunistring.so.2`, env ships `.so.5`).
Do **not** fix it by removing the package — kraken2 declares wget as a dependency and
conda uninstalls kraken2 and bracken with it. The build script shims `/usr/bin/wget`
instead.

## Databases on /media/8tb/bioThreatFinder

| Artifact | State |
|---|---|
| `taxdump/2026-09-29/` | **pinned**, md5-verified, `PROVENANCE.txt` inside. 3,016,439 nodes |
| `packs/human/refs/` | 89/89 targets, 359 genomes, 89 competitive bundles + 178 indexes, 43 GB. All bundles verified free of duplicated sequence |
| `packs/human/refs/masked/` | complete — 20.0 Mbp masked (93.9 kb vector, 2.7 kb host, 79.7 kb cross-kingdom, rest low-complexity). No genome tripped `OVER_MASKED` |
| `host/human/` | T2T-CHM13v2.0 + sr and map-ont indexes, 17 GB |
| `packs/human/markers/` | sequences intact (nt=143, aa=124, hmm=25); **manifest being rebuilt** |
| `kraken2/pluspf/` | building |
| `krakenuniq/`, `metabuli/` | not built. KrakenUniq is opt-in by design; Metabuli still pending a decision |

## The one thing that was broken at handoff

A `--only-marker` run **overwrote** `markers.manifest.tsv` with just the 19 rows it
processed, destroying the record of the 131 already built. That matters more than it
sounds: `marker_evidence.py` treats the manifest as *authoritative* for tested-vs-untested,
so every missing row silently becomes UNKNOWN — which would have taken the whole marker
tier offline while the sequence files sat intact on disk.

Fixed (subset runs now merge, and refuse to continue if the existing manifest cannot be
read) and a full pass is running to rebuild it. **On resume, confirm 163 rows.**

## Environment facts worth not rediscovering

- **rsync (port 873) is blocked** on this host; HTTPS and FTP are open. `kraken2-build`
  defaults to rsync and reports the failure as `"Error downloading assembly summary
  file"`, which looks like an NCBI outage. `--transport auto` probes and falls back.
- NCBI Entrez is being queried **anonymously** (3 req/s). Supplying an email or
  `NCBI_API_KEY` raises it to 10/s. I did not put the user's address into outbound API
  URLs without being asked.
- Nextflow self-updated to **26.04.6**, which enforces strict config/DSL syntax. Any
  pre-24 Nextflow pipeline on this machine will now fail to parse.

## Version pins

All `modules/local/*.nf` conda pins were written from memory and 11 of 13 were stale.
Now updated. Two deliberate exceptions:

- **`kraken2=2.17.2`** — not hygiene. It must match the version that *builds* the index,
  or the format may not read back.
- **`python=3.11` left alone.** Newest is 3.15.0rc2; pinning a pipeline to a release
  candidate is worse than being a version behind.

Several are major-version jumps (fastp 0.23→1.3, chopper 0.7→0.14, sylph 0.6→1.0,
mmseqs2 15→18, metabuli 1.0→1.2). **None has been run yet** — validate them with the
benchmark rather than assuming they are drop-in.

## Next steps, in order

1. Confirm the marker manifest has 163 rows, and check `bmand_18S` specifically — it is a
   `required` marker, so its failure makes *Balamuthia mandrillaris* unconfirmable, and
   sequencing is often the only route to detecting it.
2. Wait out the Kraken2 build, then verify `--report-minimizer-data` actually produces
   the distinct-minimizer and duplicity columns on a real sample.
3. Generate one small benchmark batch (`scripts/simulate_benchmark.py`) and run the
   pipeline **end to end**. That first real run will find more problems than any further
   design work — seven silent-failure bugs have been found so far and all of them shared
   the shape "produces an empty or degraded result indistinguishable from a clean
   negative".
4. Then set thresholds from measurement instead of from the literature defaults
   currently in `pack.yaml`.

## Still unresolved

- Metabuli: keep as a recall arm or drop? Independent benchmarks rank its precision
  lowest; the divergence review says it mitigates Kraken2's ~50% family-exclusion
  mis-assignment rate. Both can be true — decide with the bake-off.
- 4 near-neighbour genomes do not exist anywhere (*Giardia ardeae*, *Entamoeba coli*,
  Sandfly fever Naples phlebovirus, rat hepatitis E virus). Recorded with empty paths, so
  the weakened discrimination is visible rather than silent.
- *F. tularensis* vs *F. novicida* cannot currently be discriminated: both resolved to the
  same assembly and the duplicate was dropped. Needs a curated distinct assembly.
- `curation_required` markers (4) and `motif_todo` markers (9) need manual work; they
  cannot be automated.
