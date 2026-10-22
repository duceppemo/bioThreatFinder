# Session handoff — last updated 2026-10-21

Written because the session was ending with long builds in flight, and kept as the
running handoff since. Read this first on resume. The dated entries below are in the
order they were found, not reverse-chronological; the newest material is under
**Real data**, **Still unresolved**, and the sylph-sketch section.

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
| **Three-arm bake-off** | **DONE** 2026-10-04 → `runs/lod_v3_threearm`. Metabuli changed no finding: `findings.tsv` identical to the Kraken2-only run in every scored column, `tier` differing in 0 of 178 rows | — |
| **Divergent-agent bake-off** | **DONE** 2026-10-04 → `simulated/divergent_v1/bakeoff`. The scenarios built to favour the amino-acid arm; it changed nothing there either, and the pipeline went silent on the divergent filovirus for a reason neither arm can fix (D15) | — |
| **lod_v4 (matrix blank)** | **DONE** 2026-10-05 → `runs/lod_v4`. 24 samples, the first batch to carry a host-matrix control. `supported` 18 → 0, LoD unchanged at 1000 RPM (D14) | — |
| **D12 fine grid** | **DONE** 2026-10-06 → `runs/d12fine`. Cliff at 176.6 rpm, inside the predicted 171–184 band | — |

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
| `packs/human/refs/` | 90/90 targets, 379 manifest rows / 375 with sequence, 90 competitive bundles + 180 indexes, 43 GB. All bundles verified free of duplicated sequence. The 4 rows without sequence are near neighbours that do not exist anywhere, carried with empty paths on purpose |
| `packs/human/refs/masked/` | complete — **20,040,354 bp** over all 375 references with sequence (0.348%); vector 104,184 bp / 45 contigs, crossgroup 81,196 / 50, host 2,737 / 5. The two vacuous rows were re-masked 2026-10-04 and both shared-assembly pairs now agree and carry `SHARED_SEQUENCE` flags. Per-source byte totals in D5 were double-counted and are upper bounds, not a partition |
| `host/human/` | T2T-CHM13v2.0 + sr and map-ont indexes, 17 GB |
| `sylph/` | `human_pack_refs.syldb`, 177 MB, **375 genomes** at c=200/k=31/min-spacing=30, rebuilt 2026-10-21 with sylph 1.0.0 and a `BTF_BUILD_MANIFEST.txt`. The 2026-09-29 hand-built sketch (359 genomes) is kept beside it as `human_pack_refs.bak-2026-09-29.syldb` — the `.syldb` suffix is deliberate, `sylph inspect` refuses any other extension. Reproducible since 2026-10-05 via `scripts/build_sylph_db.sh` (6 s at 32 threads) |
| `packs/human/markers/` | 166 manifest rows, **146 built**, 33 profile HMMs (rebuilt 2026-10-21, D41). 16 can only be UNKNOWN (3 failed, 4 curation_required, 9 motif_todo) and 4 are built on fewer tiers than declared — see the vaccine/veto note below. D30 closed 11 of the 18 gaps D27 tracked |
| `kraken2/pluspf/` | **complete** 2026-10-02. `hash.k2d` 123.5 GB, `opts.k2d`, `taxo.k2d`, `seqid2taxid.map`, and `database{100,150,250}mers.kmer_distrib`. k=35, minimizer 31, full-sensitivity (no `--max-db-size`). Built against the pinned 2026-09-29 taxdump, md5-verified. Includes the 5,413 assemblies `kraken2-build` had silently dropped. `--clean` run 2026-10-03: 614 GB -> 233 GB. The remaining 125.5 GB is `database.kraken`, a Bracken BUILD intermediate - query-time Bracken reads only the `.kmer_distrib` files (measured: 0.6 s, 7.4 MB RSS). Keep it only to derive an extra read length without re-classifying the library |
| `krakenuniq/` | not built. Opt-in by design (D2) |
| `metabuli/refseq` | **built** 2026-10-03, 110 GB, NCBI-namespaced. The arm is now opt-in too (`--skip_metabuli` defaults true): measured, it changed zero findings on every scenario run, including the `divergent_*` pair meant to favour it. Kept, not removed, because its untested case — an agent absent from the index — cannot be simulated from a database. See D15 |

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

**Done 2026-10-03.** Both rows now carry real results, and the *Shigella* reference turned
out to hold **10,412 bp of vector sequence** - the Wright et al. failure mode this step
exists for:

| taxid | taxon | genome_bp | masked_bp | intervals |
|---|---|---|---|---|
| 264 | *F. tularensis* subsp. *novicida* | 1,870,206 | 24 | 1 |
| 620 | *Shigella* | 5,192,674 | 11,144 | 28 |

Each pair now agrees exactly with its assembly partner (263 and 622), which is the
invariant the bug violated, and `tests/test_pipeline.py` group `shared_assembly` pins it -
including a not-vacuous check, so a regression cannot pass by agreeing at zero.

**Fixed 2026-10-04:** re-masking both pairs in one invocation populated the flag
(`620 -> SHARED_SEQUENCE:622_Shigella_dysenteriae`, `264 -> SHARED_SEQUENCE:263_...`) and
made their per-source columns agree; `crossgroup_bp` had differed 1,835 vs 2,204 between two
rows with identical intervals. The original note follows, since the mechanism still applies
to any future subset run.

**Closed 2026-10-16 (D34), and doing it found something larger.**
`mask_pack_refs.py` used to detect sharing from RUN-TIME contig collisions, so the flag
only appeared when both partners were masked in the same invocation - the precise shape
decisions.md warns about ("a guard that asserts how a script was invoked is not a
guard"). It is now derived from the refs manifest's `accession` column, an artifact
rather than run state, and the flag is symmetric: each partner names the other.

Reading that column is what surfaced the rest. **The pack has eight shared assemblies,
not two** - six of them the Filoviridae group entry sharing an assembly with the member
taxon it was built from, which no `--only-taxid` run could ever have put in one
invocation because the partners have different taxids. **And twelve references had no
masking result at all:** every `11266_Filoviridae*` file, because the group entry was
added to the pack on 2026-10-13/14 and the masking pass ran 2026-10-03/04. Six of the
twelve had contigs in no BED anywhere, so `coverage_stats` scored them against a mask
map that did not contain them - indistinguishable from examined-and-clean.

All 371 references now have a result; the twelve came back at 0 bp masked, so the
practical cost was nil, but that is now measured rather than assumed. The script
reconciles the whole pack against the manifest on every exit path (including
`--dry-run` and the "nothing to do" early return) and exits 1 when a reference has no
result, so this cannot recur silently:

```bash
python3 scripts/mask_pack_refs.py --pack human \
        --btf-db /media/8tb/bioThreatFinder --dry-run    # read-only; exits 1 if incomplete
```

## DATA fix, closed same day: the sylph sketch predated 16 references, including all of Filoviridae

Found 2026-10-21 while reconciling the README against the built databases. `sylph/human_pack_refs.syldb`
is dated **2026-09-29** and holds **359 genomes**. The refs manifest now has 379 rows, 375 of
them with sequence, and **16 of those were downloaded after the sketch was built**:

| added | taxid | what |
|---|---|---|
| 2026-10-14 | 11266 | all **twelve** `Filoviridae` group references (D24/D25) |
| 2026-10-21 | 1639 | D42's second *L. monocytogenes* assembly |
| 2026-10-21 | 1638 / 1640 / 1642 | three *Listeria* near neighbours |

This is the familiar shape: not an error, a silent subtraction. sylph is **on by default**
and it profiles the pack, so a reference absent from the sketch is a reference that arm can
never raise — and `screen_summary.tsv` will still record the arm as `ok`, because it ran and
produced a usable report. The Filoviridae entry is the costly one: its entire purpose (D25)
is surfacing a divergent filovirus that no species-keyed row will catch, and it was measured
recovering 3 of 3 replicates **on the Kraken2 arm**. The sylph arm has never been able to
see it.

Nothing downstream is wrong, and no published finding needs revisiting — the arm was simply
blind to these taxa.

**Done 2026-10-21.** Rebuilt in 6 s at 32 threads. **sylph lives in its own env, not
`btf-refs`** — that env has no sylph at all, and the version matters because the sketch
must be readable by the version that queries it (`modules/local/detect.nf` pins
`sylph=1.0.0`, and `btf-sylph` has exactly that):

```bash
export PATH=/home/bioinfo/miniconda3/envs/btf-sylph/bin:$PATH
bash scripts/build_sylph_db.sh --pack human --btf-db /media/8tb/bioThreatFinder \
     --threads 32 --force      # --force is required; the sketch already exists
```

Verified by inspecting both sketches and diffing the genome sets, rather than by trusting
the count the builder printed:

| | old (2026-09-29) | new (2026-10-21) |
|---|---|---|
| genomes | 359 | **375** |
| `11266_Filoviridae*` | **0** | **12** |
| `1639_Listeria_monocytogenes*` | 1 | 2 |
| genomes removed | — | **none** |
| c / k / min-spacing | — | 200 / 31 / 30, matching `sylph profile`'s defaults |

The change is strictly additive: exactly the 16 expected files appeared and nothing was
dropped. The old sketch is kept as `human_pack_refs.bak-2026-09-29.syldb`; note that
`sylph inspect` **refuses a file that does not end in `.syldb`**, so a backup named
`.syldb.bak-<date>` cannot be read back and is useless as a rollback.

**The general lesson, which was D45's lesson arriving a second time:** `build_pack_refs.py`
rebuilds a competitive bundle when its inputs changed rather than when it is missing; the
sylph sketch had no such check — it is a single file that exists, so nothing noticed it was
older than the references it claimed to cover.

**Closed the same day, as D46**, in the two places that can see different things:

| where | what it compares | on disagreement |
|---|---|---|
| `main.nf` `warnIfSketchStale()`, in the sylph branch — **at launch, which is the point** | the sketch's recorded genome set against the pack, and reference mtimes against the sketch | `log.warn` naming the references, `is STALE` when the record proves it and `may be STALE` when only a timestamp suggests it |
| `scripts/build_sylph_db.sh --check` | the database's own genome list (`sylph inspect`) against `refs/`, plus mtimes, `c`/`k`/`min-spacing`, and the record against the database | names the genomes, exits **3** (0 fresh, 1 undetermined) |
| `tests/test_pipeline.py -k sketch` | the shipped sketch and its record against the shipped manifest, with no sylph needed | fails the suite |

**Why there is a record and not a cron line.** The only consumer of this sketch is a run, so
a run is where the operator is told; a schedule would report staleness on days nobody was
going to query the database. But a run cannot open a `.syldb` — the head node has no sylph,
it lives in `btf-sylph` and is used inside the process. So every build writes
`<db>.genomes.txt`, the names read back out of `sylph inspect` (never from `refs/`, which it
must be able to disagree with), and the run compares against that. `--stamp` writes it for a
database that already exists; the shipped sketch was stamped that way, 375 genomes.

Verified by pointing `--check` at the real `human_pack_refs.bak-2026-09-29.syldb`: it reports
the 16 missing references and exits 3, while the live sketch reports FRESH. The launch guard
was verified by extracting the shipped function and running it under Nextflow — against the
real manifest with that old sketch's own record, where it names all 16, and against a
synthetic pack whose database is **newer than every reference**, so mtime reports current
and only the record catches the unsketched one.

**What is open, by design:** nothing is scheduled, so a reference added and never used is
never reported — nobody launched a run. ~~The Kraken2, Metabuli and marker databases have no
freshness check at all.~~ **The marker database has one now (D47); Kraken2 and Metabuli
still do not.**

## The same predicate, third artifact: the marker database vs its catalogue (D47)

Closed 2026-10-22. The marker database's input is `packs/<pack>/markers.tsv`, so D45's
predicate transfers to it unchanged — and nothing compared them. **The shipped database was
measured first and is current** (166/166 markers, every claimed tier present, no index
older than its source), so this is a check rather than a fix.

Why it hid: `marker_evidence.py` treats `markers.manifest.tsv` as authoritative for
tested-vs-untested, so a catalogue marker with no manifest row evaluates UNKNOWN, which
never promotes a finding — the call settles a tier lower. And this pack already carries
**16 markers that can only ever be UNKNOWN** (3 failed, 4 curation_required, 9 motif_todo),
so a seventeenth is invisible among sixteen and the ready explanation is the one needing no
action.

Four comparisons, and **the mtime/set relationship is inverted from the sketch's**:

| comparison | what only it can see |
|---|---|
| `marker_id` sets | a marker added or removed with every mtime settled |
| mtime | a marker **edited** under an id the manifest already knows — D41's wrong-gene shape, which the sets cannot see |
| manifest vs combined artifact | a claimed tier the database does not hold. **The unsafe one:** downstream reports ABSENT, asserting the marker was searched for |
| index vs source | an index older than the FASTA/HMM it was built from |

```bash
# read-only; exits 0 fresh, 3 stale and names the markers, 1 could not tell.
# --db-root MATTERS: it defaults to the shared /media/36tb/db, so a bare --check reports
# this healthy pack as having no database at all.
python3 scripts/build_marker_db.py --pack human --db-root /media/8tb --check

# Or point it straight at a directory, which also checks a backup in place:
python3 scripts/build_marker_db.py --pack human --check \
        --markers-dir /media/8tb/bioThreatFinder/packs/human/markers
```

Checked at launch in `main.nf` `warnIfMarkerDbStale()` (the two TSVs only — cheap, and it
catches the operator-caused shape), in full by `--check`, and the suite runs the full check
against the shipped database on every commit. **Three things worth not rediscovering:**

- **The guard did not compile the first time.** The column parser was a closure in a local
  `def`, called as `idsOf(...)`. Nextflow 26.04 resolves `name(...)` against declared
  functions only, so `main.nf` failed to compile — and markers are ON unless
  `--skip_markers`, so that was every run. Found by extracting the function and running it
  under Nextflow, not by reading it. It is a top-level `def markerIdColumn()` now.
- **Compare presence, never sequence counts.** `concat_files()` de-duplicates within each
  marker while `n_nt_seqs` counts the per-marker archive, so they are not meant to agree.
  The count-comparing version was built and run: it calls the **correct shipped database
  stale on 116 of 166 markers**. A guard that cries stale on a good database gets switched
  off, so that version is worse than none.
- **The guard's own advice did not work, for the same reason the `--db-root` note above
  exists.** It printed `--pack human --check`, which on a healthy pack exits 3 with "no
  marker database" because of that default. The warning now passes `--markers-dir <the
  directory it judged>`, and that printed command was run verbatim to confirm it.

Verified against the retained pre-D40 database (`backups/markers-pre-D40-rebuild-20261021`,
163 rows): against the current 166-row catalogue it names exactly `salm_invA`, `salm_ttrC`,
`salm_ssaV` and exits 3, and the launch guard run under Nextflow names the same three. Suite
**364 passed, 0 failed**.

## A behaviour change to be aware of: unevaluable vaccine and veto checks now cap the tier

When this was found, cross-checking all 58 rule sets against the built manifest: 8 of 9
`vaccine_if` rules and
6 `veto_if` rules reference markers that can only ever be UNKNOWN. Under Kleene logic
UNKNOWN in those two rule types means *no downgrade* — the unsafe direction — so a Sabin
vaccine strain satisfied Enterovirus C's `confirm_requires` and reached
`confirmed_in_silico`: a vaccine strain reported as confirmed poliovirus.

An unevaluable vaccine or veto check now holds the finding at `supported`. **This reduces
confirmations**, deliberately: "supported, vaccine status untestable" costs a follow-up,
"confirmed wild-type poliovirus" costs an outbreak response. Building the 9 `motif_todo`
and 4 `curation_required` markers is what restores those confirmations.

**Re-measured 2026-10-21** against 59 rule sets and the rebuilt manifest, using
`marker_rules.compile_rule` rather than a hand grep: **9 of 10 `vaccine_if`** and **4 of
10 `veto_if`** rules still reference one of the 16 unevaluable markers. `veto_if` is down
from 6 because D30 repaired `bthai_ara`; `vaccine_if` went 8 of 9 to 9 of 10 because a
rule set was added, not because one regressed. `confirm_requires` (54 rules) and
`support_requires` (31) are **fully backed** — zero unevaluable references. The four
remaining vetoes are *F. tularensis* (`fnov_specific`), *Variola* (`vaccinia_specific`)
and Nipah/Hendra (both `henipa_cedar`).

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

1. **DONE.** The pipeline has run end to end many times (`runs/` holds eight completed
   runs). The original note follows because the ONT half is still true - no ONT data has
   ever been processed.
1b. **Run the pipeline end to end.** There is now one usable input: an ART-generated
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

## Real data — first measurement, 2026-10-08

`real/mock_v1`: ZymoBIOMICS D6300 (EVEN + LOG, Illumina and ONT) plus all ten member
organisms as single-organism Illumina libraries. 54 GB, md5-verified. Run as `zymo_rl100`,
`zymo_rl150` and `zymo_ont` because Bracken takes one read length per run.

It refuted D12's neighbour-credibility threshold (now defaulted off) and answered D14's
VERIFY on the tier ladder, cutting real-data false positives 12 → 2 with both simulated
benchmarks unchanged. See **D18**.

**`zymo_ont` COMPLETED 2026-10-16** (15h 36m, 425 CPU-hours) - the first ONT data this
project has processed end to end. `MARKER_SEARCH_AA` took **7h 20m** on one sample, so
the old 4 h default was off by nearly a factor of two and the 24 h raise was the right
size. Specificity is clean: 84 absent taxa, four findings, all held at
`attributed_to_neighbour`, zero false confirmations. Sensitivity is one of four:
*L. monocytogenes* `supported`, the other three absorbed by `attributed_to_neighbour`.

**Every ONT call is capped by the D10 identity floor and nothing else.** *Listeria* has
breadth 0.95, rpm 136,951, neighbour_ratio 9.7 and BOTH confirmatory markers PRESENT -
its only failed gate is `identity 0.8579 >= 0.9000 (ont/default)`. At 0.90 no ONT sample
can reach `confirmed_in_silico` at all. **Do not lower the floor on this evidence:** the
data is Nicholls 2019, which is R9.4.1, and D10's 0.90 was reasoned for R10.4.1 at
~97-98%. The floor was not falsified; this benchmark cannot test it. The ONT arm still
has no data on the chemistry it is calibrated for. See **D35**.

## Still unresolved

- **No real-data samplesheet declares a negative control, and the value of declaring
  one is now measured.** `runs/zymo_rl100_rescore` put two *L. monocytogenes*
  `supported` calls in the *B. subtilis* and *S. cerevisiae* isolate libraries at
  identity 0.9615/0.9646 against 0.9643 in Listeria's own library - almost certainly
  real contaminant DNA (index hopping), graded as false positives because nothing can
  tell the two apart without a control. `supported` PPV of 0.600 is therefore a FLOOR.
  Re-scoring with the yeast library declared a control (D39, 2026-10-20): PPV **0.600
  -> 0.750**, false positives **2 -> 1**, **no** per-scenario detection rate changed,
  and all 3 suppressed rows truth-absent - so the control cost nothing in sensitivity
  here. The surviving false positive clears the margin at **5.11x** against 5.0x while
  carrying Listeria's own `hly` and `prfA`: contaminating DNA, now quantitatively
  supported rather than inferred. **The decision is still open**, because the yeast
  library is another organism's library from the same run rather than a no-template
  control, and declaring it one imposes a run-wide 5x floor that on other data would
  suppress a genuine low-level positive. No samplesheet in the repo was edited;
  `benchmarks/scenarios/real_mock_communities.yaml` would need the control declared
  alongside one. See **D38**, **D39**.
  - Reaching that measurement required fixing the scorer first: `--samplesheet`
    declared controls for the QC sentence while `control_rpm` read `role` only off the
    pipeline's own rows, so declaring a control on the **re-score path** silenced
    *"NO NEGATIVE CONTROL IN RUN"*, subtracted nothing, and put *"1 negative control(s)
    in run, all clean"* on all 82 findings of a run whose declared control was itself
    reported as a `supported` *L. monocytogenes* detection at 355.5 RPM. The sheet now
    drives `role` on every row, in both directions, and a declared control absent from
    the read ledger is not credited. See **D39**.
- ~~The 0.98 illumina `min_read_identity` floor now binds on real data too.~~
  **Resolved 2026-10-21 by a second reference, not by moving the floor (D42).** The
  isolate reading 0.9643 was being compared against the only assembly the pack carried
  for it, from the other lineage. With a second assembly added: identity **0.9641 ->
  0.9947**, best-reference breadth **0.877 -> 0.982**, and **48% more reads map at all**.
  The floor stays at 0.98 and no per-taxon override exists. `--max-per-taxon` is 1 for
  bacteria and 3 for viruses by default, and the help text already gave the reason -
  intra-species diversity losing reads from a divergent strain - so this is the existing
  design extended to a bacterium on the strength of a measurement. See **D38**, **D42**.
- ~~Which other bacteria deserve a second reference is unanswered.~~ **Measured
  2026-10-21 (D43),** and without a specimen: 13,582 known assemblies across the 26
  bacterial targets, each scored against the reference the pipeline actually aligns to.
  **17 of 26 species need no addition** - one reference already reaches their most
  divergent known member - so D42's treatment is not broadly necessary. **9 do not
  span their own diversity**, and a seeded greedy cover asks for **32 further
  references**. For the D42 species itself, the single original reference would have
  failed the identity gate on **47.1% of 853 known assemblies**, falling to 3.6% with
  both: that sample was not an unlucky outlier. Nothing was added - D43 is the
  measurement and the recommendation. Two cautions recorded there: an ANI threshold
  cannot carry the verdict, because the ANI at which read identity crosses 0.98 is
  species-dependent; and one species degrades discontinuously rather than gradually, so
  its divergent members return breadth 0.002 and are non-actionable whatever the markers
  say, since the geometry gates must pass and marker rules only ever hold a finding
  down. Report: `docs/reference_span.tsv`. See **D43**.

- ~~*Salmonella* needs a real marker block, not a declaration.~~ **Done 2026-10-20
  (D40):** `salm_invA` + `salm_ttrC` + `salm_ssaV`, and *S. enterica* in the
  `ERR2935848` library is now the project's **first `confirmed_in_silico` on real data
  backed by curated markers** (markers at breadth 1.000, identity 0.985-0.989,
  475-605x; exactly one row moved of 82). invA and ttrC are both GENUS-level - S.
  bongori carries both - so the species call rests on ssaV, the SPI-2 T3SS that S.
  bongori lacks. Serovar stays out of scope per the watchlist note, and `tviB`/`viaB`
  was considered as a typhoidal proxy and rejected (also in Paratyphi C and some
  Dublin; Vi-negative Typhi exist).
- **31 of 90 watchlist taxa still have no `marker_rules.yaml` block, so
  `confirmed_in_silico` is unreachable for them by construction** (D17a's cap, which is
  correct). The binding constraint on the top tier is curation, not a threshold: the
  alignment gates are reachable on real 2x101 Illumina data (D37). The escape hatch
  `no_markers_available: true` parses (D37; it never did), but it is a statement that no
  marker exists and is the WRONG answer for a taxon with well-known targets. As of
  2026-10-20 *Orientia tsutsugamushi* is the only NOTIFIABLE BACTERIUM left without a
  block. See **D37**, **D40**.
- ~~The shipped marker database holds the wrong gene for two markers.~~ **Rebuilt
  2026-10-21 (D41):** 143 -> 146 `built`, zero status regressions, 33 profiles
  unchanged, and the three Salmonella markers now in the operational database rather
  than a scratch copy. Verified by realigning the real library against the rebuilt
  references and re-scoring: `confirmed_in_silico` fires on the production database.
  `ypest_pmt1_caf1` now holds the gene it declares. `nipah_L` no longer holds the wrong
  coding sequence - it fell back to the whole-record route, which is intended, and now
  matches its five sibling markers in strategy, coverage floor, evidence tiers and
  route. Both pending gene-tag entries were verified and removed; `GENE_TAG_MUST_FIX`
  is empty.
  - The FIRST attempt had to be rolled back: a local `hmmbuild` that resolved but could
    not run (missing shared library) failed the profile step for every marker declaring
    an HMM tier, taking `built` 143 -> 113 and `partial` 4 -> 37. Nothing to do with
    NCBI or the retrieved sequences. Restored from backup, replacement validated on one
    marker in a scratch `--db-root`, then rerun. Resolving a tool is not the same as
    running it. See **D41**.
- **16 markers that declare a gene hold a whole viral genome or segment as their
  nucleotide reference**, because the whole-record retrieval route accepts anything
  under `NT_MAX_DIRECT_LEN = 50_000`. Their `min_coverage` is 0.50, so for those rows the
  nucleotide tier asks "is half of this genome covered" rather than "is this gene
  present", and the protein and profile tiers carry the specificity. PRE-EXISTING and
  unchanged by the D41 rebuild. This is a question about what those markers should mean,
  not a bug to fix quietly, and it is recorded rather than touched. See **D41**.

- ~~Metabuli: keep as a recall arm or drop?~~ **Decided 2026-10-05: kept, defaulted OFF.**
  The bake-off ran on both the identity-matched LoD batch and the `divergent_*` pair, and
  the arm changed no finding in either — `findings.tsv` identical in every scored column,
  at ~22× Kraken2's CPU per sample. It is not removed because the one case that would
  favour it, an agent genuinely absent from the index, cannot be simulated from a
  database. Turn it on deliberately for a novel-agent question (D15).
- D12's remedy is still unchosen. Option 2 is now calibrated — neighbour "uniformity"
  (breadth ÷ uniform expectation) reads 0.114–0.124 for the absent neighbour and
  0.639–0.992 for genuinely present ones, so a threshold in the gap divides them — but no
  gate has been changed pending the fine-grid cliff measurement.
- 4 near-neighbour genomes do not exist anywhere (*Giardia ardeae*, *Entamoeba coli*,
  Sandfly fever Naples phlebovirus, rat hepatitis E virus). Recorded with empty paths, so
  the weakened discrimination is visible rather than silent.
- *F. tularensis* vs *F. novicida* cannot currently be discriminated: both resolved to the
  same assembly and the duplicate was dropped. Needs a curated distinct assembly.
- `curation_required` markers (4) and `motif_todo` markers (9) need manual work; they
  cannot be automated. These are now load-bearing: 9 of 10 `vaccine_if` rules and 4 of 10
  `veto_if` rules still reference a marker with no built tier. D30 cut the veto count from
  6 by repairing `bthai_ara`'s query (it was searching `araBAD`, an operon name NCBI does
  not index), which restores the *B. pseudomallei* and *B. mallei* vetoes. The four
  remaining retrieval gaps each carry a measured diagnosis in `markers.tsv`; note in
  particular that **`bruc_S19_eri` must not be "fixed" by retrieval** — it is a deletion,
  and supplying the eri operon would make `vaccine_if` fire on wild-type *B. abortus*.
- `-profile test` has no sample data (`tests/data/samplesheet_test.csv` is a header
  only). An empty samplesheet is now a hard error rather than a green run that analysed
  nothing, but the profile still needs a real fixture and built databases to be useful.
- ~~`screen_summary.tsv` is published but reaches neither `SCORE_EVIDENCE` nor
  `MAKE_REPORT`.~~ **Closed 2026-10-15 (D33), and it was worse than this bullet said.**
  The summary held no row at all for an arm that produced nothing, so the file whose
  job was to be the trace of an ignored failure recorded that failure as absence -
  indistinguishable from a *disabled* arm, because `fillMissing()` pads both to `[]`.
  The workflow now passes `--expected-arms`, the four states are named (`ok`,
  `empty_report`, `no_output`, `not_run`), and the summary reaches `MAKE_REPORT`,
  which banners a dead arm. It still does NOT reach `SCORE_EVIDENCE`, deliberately:
  it is coverage, and must not move a call. The same fix surfaced the unclassified
  fraction, which had been computed, logged, written and read by nothing - and made
  `--high_unclassified` reachable, which it had never been. **Corrected 2026-10-19
  (D36):** the five states became five, not four - `not_applicable` joined them on the
  first ONT run (D34/e38372e) - and `empty_report` was raising a false alarm. sylph
  searches the PACK and answers a sample with no watchlist organism in it with an
  empty profile, so four of eleven correct answers in the first real Illumina
  re-score were bannered as a classifier that had produced nothing. An empty report is
  now an alarm only for the arms that search a comprehensive database; the run-level
  check that replaces it - an arm empty on EVERY sample, which is what a truncated
  database looks like - lives in `make_report.py`, the only place that sees every
  sample at once.
