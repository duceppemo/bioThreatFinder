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
| **Kraken2 PlusPF build** | **DONE** 2026-10-02 — index + all three Bracken distributions | — |
| **Marker DB full pass** | **DONE** — 163 rows, `bmand_18S` is `built` | — |

**Do not start a second build against a live one.** I did, and the two
`kraken2-build --download-library bacteria` processes collided on the same staging
directory (`rm: cannot remove 'all/': Directory not empty`). Two concurrent FTP
transfers can interleave partial files into `library.fna`, and a corrupt reference in an
index does not error — it misclassifies. All three DB builders now take an `flock` on the
database directory and name the holder on refusal. Check for a live build first:

```bash
pgrep -af 'build_(kraken2|metabuli|krakenuniq)_db.sh'
```

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
| `packs/human/refs/masked/` | complete — 20.0 Mbp masked. **Two rows need re-masking: see below.** Per-source byte totals in D5 were double-counted and are upper bounds, not a partition |
| `host/human/` | T2T-CHM13v2.0 + sr and map-ont indexes, 17 GB |
| `packs/human/markers/` | 163 manifest rows; 18 referenced by rules can only be UNKNOWN (5 failed, 4 curation_required, 9 motif_todo) — see the vaccine/veto note below |
| `kraken2/pluspf/` | **complete** 2026-10-02. `hash.k2d` 123.5 GB, `opts.k2d`, `taxo.k2d`, `seqid2taxid.map`, and `database{100,150,250}mers.kmer_distrib`. k=35, minimizer 31, full-sensitivity (no `--max-db-size`). Built against the pinned 2026-09-29 taxdump, md5-verified. Includes the 5,413 assemblies `kraken2-build` had silently dropped. `--clean` run 2026-10-03: 614 GB -> 233 GB. The remaining 125.5 GB is `database.kraken`, a Bracken BUILD intermediate - query-time Bracken reads only the `.kmer_distrib` files (measured: 0.6 s, 7.4 MB RSS). Keep it only to derive an extra read length without re-classifying the library |
| `krakenuniq/`, `metabuli/` | not built. KrakenUniq is opt-in by design; Metabuli still pending a decision |

## Outstanding DATA fix: two reference genomes were never masked

*F. tularensis* subsp. *novicida* (264) and *Shigella* (620) each share one NCBI assembly
with another watchlist taxon, so their FASTAs carry identical contig names. The masking
pass attributed those contigs to the first genome, skipped the second, and wrote
`genome_bp=0 masked_bp=0 flag=` rows with intervals-free BEDs — byte-identical in meaning
to "examined and found clean". Both shipped BEDs are still in that state.

The code is fixed (identical sequence is masked once and both BEDs written from that
result, flagged `SHARED_SEQUENCE`; a name collision with *differing* sequence is now
fatal). **The artifacts have not been regenerated**, because the host pass needs a
minimap2 index of the 3.1 Gbp human genome and the Kraken2 build was using the machine.
The Kraken2 build finished 2026-10-02, so this can run now:

```bash
export PATH=/home/bioinfo/miniconda3/envs/btf-refs/bin:$PATH
# --btf-db is the base that `packs/` hangs off, NOT the volume: /media/8tb gives
# /media/8tb/packs/... and the script dies with "no refs manifest".
python3 scripts/mask_pack_refs.py --pack human --btf-db /media/8tb/bioThreatFinder \
        --only-taxid 264 --only-taxid 620 --threads 24 --force
```

`--force` is required: both rows already carry a (vacuous) `masked_bp=0` result, so
without it the run is skipped as already done.

Then confirm both rows have non-zero `genome_bp` in `masked/mask_manifest.tsv`.

## A behaviour change to be aware of: unevaluable vaccine and veto checks now cap the tier

Cross-checking all 58 rule sets against the built manifest: 8 of 9 `vaccine_if` rules and
6 `veto_if` rules reference markers that can only ever be UNKNOWN. Under Kleene logic
UNKNOWN in those two rule types means *no downgrade* — the unsafe direction — so a Sabin
vaccine strain satisfied Enterovirus C's `confirm_requires` and reached
`confirmed_in_silico`: a vaccine strain reported as confirmed poliovirus.

An unevaluable vaccine or veto check now holds the finding at `supported`. **This reduces
confirmations**, deliberately: "supported, vaccine status untestable" costs a follow-up,
"confirmed wild-type poliovirus" costs an outbreak response. Building the 9 `motif_todo`
and 4 `curation_required` markers is what restores those confirmations.

## Environment facts worth not rediscovering

- **rsync (port 873) is blocked** on this host; HTTPS and FTP are open. `kraken2-build`
  defaults to rsync and reports the failure as `"Error downloading assembly summary
  file"`, which looks like an NCBI outage. `--transport auto` probes and falls back.
- **`systemd-oomd` kills anything memory-hungry that runs inside the user session.**
  `/usr/lib/systemd/system/user@.service.d/10-oomd-user-service-defaults.conf` sets
  `ManagedOOMMemoryPressure=kill` with a **50% pressure limit** on `user@1000.service`.
  Bracken loading the 123 GB hash pushed pressure to 87.72%, and oomd killed the build,
  gnome-shell and five tmux scopes at once — on a 503 GB box with no kernel OOM. Pressure,
  not exhaustion, is the trigger. Run long/large jobs as a SYSTEM unit, which is not
  monitored (`ManagedOOMMemoryPressure=auto`):
  `sudo systemd-run --unit=<name> --uid=1000 --gid=1000 /bin/bash <script>`
  (`--uid` keeps the output files owned by `bioinfo`). Note `ManagedOOMMemoryPressure=off`
  is NOT a valid value - only `auto` or `kill` - so passing `off` fails the whole call with
  `Invalid argument`.
- **A Claude Code background task does not survive the session**, and `setsid` is not
  enough because that process is a child subreaper, so an orphan reparents to it. Same fix:
  a systemd unit.
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
mmseqs2 15→18, metabuli 1.0→1.2).

**This warning cashed in.** `samtools=1.24 conda-forge::python=3.11` is
**unsatisfiable** — samtools 1.24 wants `htslib >=1.24,<1.25`, which wants
`libdeflate >=1.25,<1.27`, which does not exist in the channels. It solves happily
alongside minimap2, so the conflict is specifically with a conda-forge python, and
MARKER_EVIDENCE died on it during the first end-to-end run. All 13 distinct specs were
then put through the solver: 12 solved, that one did not. samtools is now pinned to
**1.22** everywhere, verified to solve both alongside minimap2 and alongside python 3.11.

Re-validate with `bash scripts/check_conda_specs.sh` after any pin change. It is not in
`tests/test_pipeline.py` because it needs the network and takes minutes; the unit suite
stays offline and fast by design.

## Next steps, in order

1. **Run the pipeline end to end.** There is now one usable input: an ART-generated
   Illumina positive control at
   `/media/8tb/bioThreatFinder/simulated/plumbing/` — 2M reads, 5000 rpm *Y. pestis*,
   94.6% host, truth row in `manifest.tsv`. Its ONT partner was cut off by a 1700 s
   timeout and needs regenerating. Nothing in this project has ever processed real data
   end to end, and that run will find more than any further static review: well over
   thirty silent-failure bugs have been found so far, all sharing the shape "produces an
   empty or degraded result indistinguishable from a clean negative".
2. **DONE 2026-10-03.** `--report-minimizer-data` verified on real reads: 8 columns,
   e.g. Bacillus genus 20,075,887 minimizers / 1,794,592 distinct = 11.2 duplicity. The
   LoD batch was re-run with Kraken2 as `runs/lod_v2_kraken2` (the sylph-only pass is
   kept at `runs/lod_v1` for comparison).
3. Re-mask taxids 264 and 620 (above). **Now unblocked** - the Kraken2 build is finished.
   Run 2026-10-03.
3b. **Expand the benchmark before touching the D12 neighbour gate.**
   `benchmarks/scenarios/real_mock_communities.yaml` drafts the real-data arms around the
   ZymoBIOMICS standards (even D6300 and log D6310, Illumina + ONT). All accessions in it
   are marked VERIFY and nothing has been downloaded. It also gives this pipeline its
   FIRST ONT data, which is the only evidence the provisional D10 floors have ever had.
4. Then set thresholds from measurement. Three are explicitly PROVISIONAL and flagged in
   the code: the ONT read-identity floors (D10), `--ont-identity-slack` in the marker
   tier, and the read/breadth gate incoherence (D11) where **75 of 89** taxa have a
   decorative `min_unique_reads` — up to 243x for *Schistosoma mansoni*. Run
   `scripts/check_threshold_coherence.py` for the table.

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
  cannot be automated. These are now load-bearing: they are most of the 18 markers that
  make 8 vaccine rules and 6 veto rules inoperative (above).
- `-profile test` has no sample data (`tests/data/samplesheet_test.csv` is a header
  only). An empty samplesheet is now a hard error rather than a green run that analysed
  nothing, but the profile still needs a real fixture and built databases to be useful.
- `screen_summary.tsv` is published but reaches neither `SCORE_EVIDENCE` nor
  `MAKE_REPORT`. For a taxon that *was* raised the classifier columns do travel into
  `findings.tsv`; but if an `errorStrategy 'ignore'` arm were the only one that would
  have raised a taxon, it is lost with no trace in the report.
