# Operations — last updated 2026-10-22

How to build, rebuild and verify the databases this pipeline needs, and the environment
traps that cost a day each to find. This file was `docs/STATUS.md`, a session handoff, until
2026-10-22 (D50); its content is operational knowledge that was paid for once and should not
be paid for twice, so it was reframed rather than trimmed.

**Two kinds of material, labelled throughout.** Unlabelled sections are portable: true for
any deployment. Sections marked **HOST-SPECIFIC** record what was observed on the
development machine — a 503 GB Ubuntu workstation with one fast NVMe volume and one large
slow one — and may not hold on yours. Nothing in a host-specific section is configuration
the pipeline reads, and no path from that machine survives anywhere in the repo outside the
dated entries in `docs/decisions.md` (D51).

**Two roots, set once.** The pipeline requires both and refuses to start without them,
naming both: `$BTF_DB_ROOT` (`--db_root`) for large shared read-only resources, and
`$BTF_DB` (`--btf_db`) for project-owned artifacts on the fastest disk available.

```bash
export BTF_DB_ROOT=/path/to/shared
export BTF_DB=/path/to/project/bioThreatFinder     # the last segment is required
```

**Every build script and every check below derives its own paths from `$BTF_DB`, and
refuses to run when it is unset** (D51) — they used to default to the two volumes above,
so a clone built against paths it did not have. Two flags still override, and they do not
mean the same thing:

| In a build script | Means | Example, for `BTF_DB=/fast/bioThreatFinder` |
|---|---|---|
| `--btf-db` | exactly `$BTF_DB` — the base that `packs/` hangs off, **not** the volume. Passing the volume gives `<volume>/packs/...` and the script dies with "no refs manifest" | `/fast/bioThreatFinder` |
| `--db-root` | the **parent** of `$BTF_DB`: these scripts append `bioThreatFinder/` themselves | `/fast` |

`$BTF_DB` must end in `bioThreatFinder` for that second form to land where the pipeline
reads. Neither flag is needed once the variable is exported, and **`--db-root` is not
`$BTF_DB_ROOT`** — passing the shared root to it is the mistake it used to default to.
Environments are entered with `conda activate <name>` from `envs/*.yml`.

## Read-only checks, all of them

None of these write anything. The three database checks share an exit convention — **0**
the answer is yes, **3** no and it names what is wrong, **1** it could not tell — which
distinguishes a bad answer from a failure to get one. Every command and exit code in this
table was run against the shipped databases on 2026-10-22.

| Check | Command | Exit | What it answers |
|---|---|---|---|
| Competitive bundles (D45) | applied by `build_pack_refs.py` itself when it decides what to rebuild | — | is each bundle's recorded input set, and each input's mtime, still what it was built from |
| sylph sketch (D46) | `bash scripts/build_sylph_db.sh --pack human --check` | 0/3/1 | does the sketch cover every pack reference |
| Marker database (D47) | `python3 scripts/build_marker_db.py --pack human --check` | 0/3/1 | does the database cover `markers.tsv`, on every tier it claims |
| Kraken2 index (D48) | `python3 scripts/check_index_coverage.py --pack human` | 0/3/1 — **3 today**, by design, see D48 | can the index report every watchlist taxon **at all** |
| Masking completeness (D34) | `python3 scripts/mask_pack_refs.py --pack human --dry-run` | 0, or 1 if a reference has no result | does every reference have a masking result |
| Threshold coherence | `python3 scripts/check_threshold_coherence.py` | report, always 0 | which per-taxon gates can never bind |
| Conda specs | `bash scripts/check_conda_specs.sh` | 0 or nonzero | does every per-process pin still solve (network, minutes) |
| The suite | `python3 tests/test_pipeline.py` | 0 or nonzero | 391 checks with samtools on `PATH` and `$BTF_DB` set; 330 pass and 14 skip on a bare clone with neither. No network, ~30 s |

**Why there are this many.** Well over thirty bugs found in this project so far share one
shape: *produces an empty or degraded result indistinguishable from a clean negative.* A
database that opens, an arm that runs, a report that parses — none of those distinguish a
correct answer from a silent subtraction. Each check above exists because something in that
list got past a check that only asked whether the artifact was there.

Only the sylph check needs a tool: `conda activate btf-sylph` first, or it exits **1**
("sylph not found"), which is the undetermined code and not a stale database. The rest are
stdlib Python and run from a bare clone — deliberately, since the head node has neither
sylph nor kraken2 and the checks that matter at launch must work without them.

**Two gaps, both known.** Metabuli has no freshness or coverage check at all; it is opt-in
and off by default (`--skip_metabuli` defaults true), which is why it has stayed unchecked.
And the bundle predicate (D45) is the one with no standalone read-only form — it is a
rebuild decision inside the builder, so asking "are my bundles current?" means invoking the
builder, whose `--dry-run` reaches the network. The answer is still correct whenever a build
runs, which is where a repair would happen anyway.

**Three of these also run at launch, which is the point:** the only consumer of any of these
databases is a run, so a run is where the operator is told. A schedule would report staleness
on days nobody was going to query the database. `main.nf` carries `warnIfSketchStale()`,
`warnIfMarkerDbStale()` and `warnIfPackPartlyMasked()` — each a `log.warn`, never a refusal,
because a stale database still answers and the operator may know why. **Kraken2 index
coverage is on demand only**, deliberately: its answer is `3` today and will stay `3` until
someone adds a helminth library (D48), so warning at every launch would train the operator
to ignore the warnings that matter.

## Build provenance, and how to resume an interrupted build

Every artifact below is built. The dates are what the shipped databases were built from;
a rebuild should reproduce them, and a figure that disagrees is the thing to investigate.

| Build | Outcome |
|---|---|
| **Kraken2 PlusPF** | 2026-10-02 — index plus all three Bracken distributions |
| **Marker DB full pass** | 163 rows, `bmand_18S` is `built`; superseded by the D41 rebuild to 166 rows |
| **Three-arm bake-off** | 2026-10-04 → `runs/lod_v3_threearm`. Metabuli changed no finding: `findings.tsv` identical to the Kraken2-only run in every scored column, `tier` differing in 0 of 178 rows |
| **Divergent-agent bake-off** | 2026-10-04 → `simulated/divergent_v1/bakeoff`. The scenarios were built to favour the amino-acid arm; it changed nothing there either, and the pipeline went silent on the divergent filovirus for a reason neither arm can fix (D15) |
| **lod_v4 (matrix blank)** | 2026-10-05 → `runs/lod_v4`. 24 samples, the first batch to carry a host-matrix control. `supported` 18 → 0, LoD unchanged at 1000 RPM (D14) |
| **D12 fine grid** | 2026-10-06 → `runs/d12fine`. Cliff at 176.6 rpm, inside the predicted 171–184 band |

**Never start a second build against a live one.** This was done once: two
`kraken2-build --download-library bacteria` processes collided on the same staging
directory (`rm: cannot remove 'all/': Directory not empty`). Two concurrent FTP transfers
can interleave partial files into `library.fna`, and **a corrupt reference in an index does
not error — it misclassifies.** All three DB builders now take an `flock` on the database
directory and name the holder on refusal. Check for a live build first:

```bash
pgrep -af 'build_(kraken2|metabuli|krakenuniq)_db.sh'
```

Both the Kraken2 and the marker build are resumable and skip completed work:

```bash
conda activate btf-kraken2
bash scripts/build_kraken2_db.sh --taxdump "$BTF_DB/taxdump/2026-09-29" --threads 40

conda activate btf-refs
python3 scripts/build_marker_db.py --pack human --taxdump "$BTF_DB/taxdump/2026-09-29" \
     --threads 24
```

A long build should not be a child of an interactive session — see **Environment facts**
below for why `setsid` and background tasks are not enough.

## Environments — three of them, and the split is not tidiness

| Env | Holds | Note |
|---|---|---|
| `btf-refs` | datasets, minimap2, seqkit, samtools, mmseqs2, hmmer, mafft | the deployed copy on the development host also carries kraken2 2.0.8-beta from the mistake below — do not use it for kraken2. `envs/btf-refs.yml` deliberately does not list kraken2 at all |
| `btf-kraken2` | **kraken2 2.17.2**, bracken | deliberately separate |
| `btf-sylph` | **sylph 1.0.0** | matches the pin in `modules/local/detect.nf` |

**Why separate:** installing kraken2 into `btf-refs` alongside mmseqs2/hmmer/mafft made the
solver resolve **kraken2 2.0.8-beta (2019)**, which has no `--report-minimizer-data` — the
flag the whole stage-1 design rests on. An index built with it would work and carry none of
the evidence the pipeline scores on. The rest of `btf-refs` was audited afterwards:
everything else was current, so this was specific to kraken2's perl/wget dependencies.
`build_kraken2_db.sh` now refuses anything below 2.1.0 and verifies the flag is actually
advertised. This is also why `envs/*.yml` are hand-written rather than
`conda env export` output: an export of the deployed `btf-refs` would ship the trap as the
recommended install.

**Why `btf-sylph` is pinned, not merely separate:** a sketch must be readable by the
version that queries it, and `modules/local/detect.nf` pins `sylph=1.0.0`, so a sketch
built by another major version is a database the run cannot open. The 0.6 → 1.0 jump was a
major version change. It is also why the sketch's freshness check (D46) reads a sidecar
record rather than the database: sylph lives only in that environment, so the Nextflow head
node has no `sylph inspect`.

**HOST-SPECIFIC:** `btf-refs`' wget is broken (linked against `libunistring.so.2`, the env
ships `.so.5`). Do **not** fix it by removing the package — kraken2 declares wget as a
dependency and conda uninstalls kraken2 and bracken with it. The build script shims
`/usr/bin/wget` instead.

## Database roots: what refuses, and what that buys

Every script that writes or reads a built database derives its paths from `$BTF_DB` and
**refuses to run when it cannot** (D51). The refusal names the variable, the flag and a
runnable example, because the failure it replaces was not a crash:

| Was | Is |
|---|---|
| a build wrote to a fixed path on a machine that has no such volume | refusal naming `$BTF_DB`, before anything is downloaded |
| `--check` read a root that was not yours and reported a **healthy pack as having no database** | the default is the pack you built, so the bare form is correct |
| an empty root silently joined into `/bioThreatFinder/packs/human/refs` — a real-looking path that exists nowhere | refused, with the two ways to answer |

Three checks take the artifact's path directly and so work with **no environment at all**:
`build_marker_db.py --markers-dir`, `check_index_coverage.py --kraken2-db`, and
`build_sylph_db.sh --refs-dir` + `--out`. That is deliberate: a check that needs a
configured deployment is a check nobody runs on a fresh clone, and these are the ones a
launch warning prints for someone to paste.

`tests/test_pipeline.py -k roots` pins both halves — no host path in any script, and a
refusal that names `$BTF_DB` in each of the eight builders. The first version of that scan
only looked at assignments, and a mutation test walked through it: a default reached
through a conditional expression is still a default.

## The built databases — reference state

Inventory under `$BTF_DB`, as built on the development host. The figures are the reference
state: this is what a correct rebuild produces, so a figure that disagrees is the thing to
investigate.

| Artifact | State |
|---|---|
| `taxdump/2026-09-29/` | **pinned**, md5-verified, `PROVENANCE.txt` inside. 3,016,439 nodes |
| `packs/human/refs/` | 90/90 targets, 379 manifest rows / 375 with sequence, 90 competitive bundles + 180 indexes, 43 GB. All bundles verified free of duplicated sequence. The 4 rows without sequence are near neighbours that do not exist anywhere, carried with empty paths on purpose |
| `packs/human/refs/masked/` | complete — **20,040,354 bp** over all 375 references with sequence (0.348%); vector 104,184 bp / 45 contigs, crossgroup 81,196 / 50, host 2,737 / 5. The two vacuous rows were re-masked 2026-10-04 and both shared-assembly pairs now agree and carry `SHARED_SEQUENCE` flags. Per-source byte totals in D5 were double-counted and are upper bounds, not a partition |
| `host/human/` | T2T-CHM13v2.0 + sr and map-ont indexes, 17 GB |
| `sylph/` | `human_pack_refs.syldb`, 177 MB, **375 genomes** at c=200/k=31/min-spacing=30, rebuilt 2026-10-21 with sylph 1.0.0 and a `BTF_BUILD_MANIFEST.txt`. The 2026-09-29 hand-built sketch (359 genomes) is kept beside it as `human_pack_refs.bak-2026-09-29.syldb` — the `.syldb` suffix is deliberate, `sylph inspect` refuses any other extension. Reproducible since 2026-10-05 via `scripts/build_sylph_db.sh` (6 s at 32 threads) |
| `packs/human/markers/` | 166 manifest rows, **146 built**, 33 profile HMMs (rebuilt 2026-10-21, D41). 16 can only be UNKNOWN (3 failed, 4 curation_required, 9 motif_todo) and 4 are built on fewer tiers than declared — see the vaccine/veto note below |
| `kraken2/pluspf/` | **complete** 2026-10-02. `hash.k2d` 123.5 GB, `opts.k2d`, `taxo.k2d`, `seqid2taxid.map`, and `database{100,150,250}mers.kmer_distrib`. k=35, minimizer 31, full-sensitivity (no `--max-db-size`). Built against the pinned 2026-09-29 taxdump, md5-verified. Includes the 5,413 assemblies `kraken2-build` had silently dropped. `--clean` run 2026-10-03: 614 GB → 233 GB. The remaining 125.5 GB is `database.kraken`, a Bracken BUILD intermediate — query-time Bracken reads only the `.kmer_distrib` files (measured: 0.6 s, 7.4 MB RSS). Keep it only to derive an extra read length without re-classifying the library |
| `krakenuniq/` | not built. Opt-in by design (D2) |
| `metabuli/refseq` | **built** 2026-10-03, 110 GB, NCBI-namespaced. The arm is opt-in too (`--skip_metabuli` defaults true): measured, it changed zero findings on every scenario run, including the `divergent_*` pair meant to favour it. Kept, not removed, because its untested case — an agent absent from the index — cannot be simulated from a database. See D15 |

**What is re-creatable and what is not.** The refs (43 GB) and markers (61 MB) are
re-fetchable: both manifests pin accession **and** md5, so a rebuild is verifiable rather
than merely plausible. The Kraken2 (233 GB) and Metabuli (110 GB) indexes are
reproducible from the scripts but impractical to move. **The pinned taxdump is the
exception:** NCBI overwrites `new_taxdump.tar.gz` in place, so that exact snapshot cannot
be re-downloaded, and the watchlist resolution, the Kraken2 taxonomy and every taxid in
`findings.tsv` key to it.

## Data fix, closed: two reference genomes were never masked

*F. tularensis* subsp. *novicida* (264) and *Shigella* (620) each share one NCBI assembly
with another watchlist taxon, so their FASTAs carry identical contig names. The masking
pass attributed those contigs to the first genome, skipped the second, and wrote
`genome_bp=0 masked_bp=0 flag=` rows with intervals-free BEDs — byte-identical in meaning
to "examined and found clean".

The code is fixed: identical sequence is masked once and both BEDs written from that
result, flagged `SHARED_SEQUENCE`; a name collision with *differing* sequence is now fatal.
Regenerating the artifacts was blocked on the host pass needing a minimap2 index of the
3.1 Gbp human genome while the Kraken2 build held the machine, and ran once that finished:

```bash
conda activate btf-refs
python3 scripts/mask_pack_refs.py --pack human --btf-db "$BTF_DB" \
        --only-taxid 264 --only-taxid 620 --threads 24 --force
```

`--force` is required: both rows already carry a (vacuous) `masked_bp=0` result, so without
it the run is skipped as already done.

**Done 2026-10-03.** Both rows now carry real results, and the *Shigella* reference turned
out to hold **10,412 bp of vector sequence** — the Wright et al. failure mode this step
exists for:

| taxid | taxon | genome_bp | masked_bp | intervals |
|---|---|---|---|---|
| 264 | *F. tularensis* subsp. *novicida* | 1,870,206 | 24 | 1 |
| 620 | *Shigella* | 5,192,674 | 11,144 | 28 |

Each pair now agrees exactly with its assembly partner (263 and 622), which is the
invariant the bug violated, and `tests/test_pipeline.py` group `shared_assembly` pins it —
including a not-vacuous check, so a regression cannot pass by agreeing at zero.

**Fixed 2026-10-04:** re-masking both pairs in one invocation populated the flag
(`620 -> SHARED_SEQUENCE:622_Shigella_dysenteriae`, `264 -> SHARED_SEQUENCE:263_...`) and
made their per-source columns agree; `crossgroup_bp` had differed 1,835 vs 2,204 between two
rows with identical intervals. The note above is kept because the mechanism still applies to
any future subset run.

**Closed 2026-10-16 (D34), and doing it found something larger.**
`mask_pack_refs.py` used to detect sharing from RUN-TIME contig collisions, so the flag only
appeared when both partners were masked in the same invocation — the precise shape
decisions.md warns about ("a guard that asserts how a script was invoked is not a guard").
It is now derived from the refs manifest's `accession` column, an artifact rather than run
state, and the flag is symmetric: each partner names the other.

Reading that column is what surfaced the rest. **The pack has eight shared assemblies, not
two** — six of them the Filoviridae group entry sharing an assembly with the member taxon it
was built from, which no `--only-taxid` run could ever have put in one invocation because
the partners have different taxids. **And twelve references had no masking result at all:**
every `11266_Filoviridae*` file, because the group entry was added to the pack on
2026-10-13/14 and the masking pass ran 2026-10-03/04. Six of the twelve had contigs in no
BED anywhere, so `coverage_stats` scored them against a mask map that did not contain them —
indistinguishable from examined-and-clean.

All 371 references now have a result; the twelve came back at 0 bp masked, so the practical
cost was nil, but that is now measured rather than assumed. The script reconciles the whole
pack against the manifest on every exit path (including `--dry-run` and the "nothing to do"
early return) and exits 1 when a reference has no result, so this cannot recur silently.

## Data fix, closed same day: the sylph sketch predated 16 references, including all of Filoviridae

Found 2026-10-21 while reconciling the README against the built databases.
`sylph/human_pack_refs.syldb` was dated **2026-09-29** and held **359 genomes**. The refs
manifest has 379 rows, 375 of them with sequence, and **16 of those were downloaded after
the sketch was built**:

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
recovering 3 of 3 replicates **on the Kraken2 arm**. The sylph arm had never been able to
see it. Nothing downstream was wrong, and no published finding needed revisiting — the arm
was simply blind to these taxa.

**Done 2026-10-21.** Rebuilt in 6 s at 32 threads:

```bash
conda activate btf-sylph
bash scripts/build_sylph_db.sh --pack human --btf-db "$BTF_DB" \
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

**Closed the same day, as D46**, in the three places that can see different things:

| where | what it compares | on disagreement |
|---|---|---|
| `main.nf` `warnIfSketchStale()`, in the sylph branch — **at launch, which is the point** | the sketch's recorded genome set against the pack, and reference mtimes against the sketch | `log.warn` naming the references, `is STALE` when the record proves it and `may be STALE` when only a timestamp suggests it |
| `scripts/build_sylph_db.sh --check` | the database's own genome list (`sylph inspect`) against `refs/`, plus mtimes, `c`/`k`/`min-spacing`, and the record against the database | names the genomes, exits **3** |
| `tests/test_pipeline.py -k sketch` | the shipped sketch and its record against the shipped manifest, with no sylph needed | fails the suite |

**Why there is a record and not a cron line.** The only consumer of this sketch is a run, so
a run is where the operator is told; a schedule would report staleness on days nobody was
going to query the database. But a run cannot open a `.syldb` — the head node has no sylph,
it lives in `btf-sylph` and is used inside the process. So every build writes
`<db>.genomes.txt`, the names read back out of `sylph inspect` (never from `refs/`, which it
must be able to disagree with), and the run compares against that. `--stamp` writes it for a
database that already exists; the shipped sketch was stamped that way, 375 genomes.

Verified by pointing `--check` at the real `human_pack_refs.bak-2026-09-29.syldb`: it
reports the 16 missing references and exits 3, while the live sketch reports FRESH. The
launch guard was verified by extracting the shipped function and running it under Nextflow —
against the real manifest with that old sketch's own record, where it names all 16, and
against a synthetic pack whose database is **newer than every reference**, so mtime reports
current and only the record catches the unsketched one.

**What is open, by design:** nothing is scheduled, so a reference added and never used is
never reported — nobody launched a run. The marker database has its own check (D47).
Kraken2's question turned out to be coverage rather than freshness and is answered in D48.
Metabuli is still unchecked, and is off by default.

## The Kraken2 index: 19 of 90 watchlist taxa have no node in it (D48)

Found 2026-10-22 by asking the question D47 said was the right one for this index, whose
inputs are NCBI library snapshots rather than `packs/`: **not "is it stale" but "can it
report every watchlist taxon at all".** A taxon with no node in the pruned taxonomy can
never appear in a Kraken2 report — the arm runs, `screen_summary.tsv` records it `ok`, and
the row simply never exists.

```bash
python3 scripts/check_index_coverage.py --pack human
```

**19 of 90 are absent:** the helminths, several protozoa, and most of the dimorphic and
opportunistic fungi — 16 notifiable, 2 emerging, 1 high_consequence. This is **not
staleness and not a build error**: it follows from PlusPF's library selection (no helminth
library exists in it at all), so *rebuilding against the same libraries changes nothing*.
Do not spend a day and 123 GB on that.

**Why it is tolerable, and the part that must stay true:** every one of the 19 has a pack
reference and is in the sylph sketch, and sylph is on by default with
`combining_rule: any_classifier_or`. So stage 1 still reaches all 19 — **one arm is blind,
not the stage.** A taxon with no Kraken2 node *and* no sketched reference would be invisible
to stage 1 entirely; the check reports that separately and the suite fails if it appears.

It also sharpens D46 in hindsight: for these 19, sylph is not a redundant second opinion, it
is the only stage-1 arm that can see them — so the 22 days the sketch was stale cost more for
this group than for any other.

**`taxo.k2d` is read with no kraken2** (the head node has none — D46's constraint again):
`K2TAXDAT` magic, three uint64 counts, then 7×uint64 per node with the taxid at field 5.
`32 + nodes*56 + name_len + rank_len` must equal the file size **to the byte** or the reader
refuses — a wrong stride does not fail, it reads the wrong 8 bytes of every record and
reports them as taxids confidently. Two overclaims are blocked by tests: the pruned taxonomy
(62,771 nodes) is never compared to the pinned taxdump's 3,016,439, and a present node is
called necessary, not sufficient.

**The open decision:** `pack.yaml`'s `eukaryote_strategy: kraken2_then_align` nominates an
arm that does not contain most of the eukaryotes. The comment that claimed otherwise has
been corrected; the value is unchanged because nothing reads the key. Adding a helminth
library has to be weighed against stating sylph + stage-2 alignment as the intended
eukaryote route. See D48.

## The same predicate, third artifact: the marker database vs its catalogue (D47)

Closed 2026-10-22. The marker database's input is `packs/<pack>/markers.tsv`, so D45's
predicate transfers to it unchanged — and nothing compared them. **The shipped database was
measured first and is current** (166/166 markers, every claimed tier present, no index older
than its source), so this is a check rather than a fix.

Why it hid: `marker_evidence.py` treats `markers.manifest.tsv` as authoritative for
tested-vs-untested, so a catalogue marker with no manifest row evaluates UNKNOWN, which
never promotes a finding — the call settles a tier lower. And this pack already carries **16
markers that can only ever be UNKNOWN** (3 failed, 4 curation_required, 9 motif_todo), so a
seventeenth is invisible among sixteen and the ready explanation is the one needing no
action.

Four comparisons, and **the mtime/set relationship is inverted from the sketch's**:

| comparison | what only it can see |
|---|---|
| `marker_id` sets | a marker added or removed with every mtime settled |
| mtime | a marker **edited** under an id the manifest already knows — D41's wrong-gene shape, which the sets cannot see |
| manifest vs combined artifact | a claimed tier the database does not hold. **The unsafe one:** downstream reports ABSENT, asserting the marker was searched for |
| index vs source | an index older than the FASTA/HMM it was built from |

```bash
python3 scripts/build_marker_db.py --pack human --check

# Or point it straight at a directory, which also checks a backup in place and needs no
# environment at all:
python3 scripts/build_marker_db.py --pack human --check \
        --markers-dir "$BTF_DB/packs/human/markers"
```

**`--db-root` used to matter here and no longer does.** It defaulted to the shared root,
where this pack has no marker database, so a bare `--check` reported a healthy pack as
having none — exit 3 from the command the launch warning had just printed. Its default now
derives from `$BTF_DB` (D51), so the bare form is correct; the warning still prints
`--markers-dir`, because that works with no environment at all.

Checked at launch in `main.nf` `warnIfMarkerDbStale()` (the two TSVs only — cheap, and it
catches the operator-caused shape), in full by `--check`, and the suite runs the full check
against the shipped database on every commit. **Three things worth not rediscovering:**

- **The guard did not compile the first time.** The column parser was a closure in a local
  `def`, called as `idsOf(...)`. Nextflow 26.04 resolves `name(...)` against declared
  functions only, so `main.nf` failed to compile — and markers are ON unless
  `--skip_markers`, so that was every run. Found by extracting the function and running it
  under Nextflow, not by reading it. It is a top-level `def markerIdColumn()` now, and CI
  compiles `main.nf` so that reading is never the only check again.
- **Compare presence, never sequence counts.** `concat_files()` de-duplicates within each
  marker while `n_nt_seqs` counts the per-marker archive, so they are not meant to agree.
  The count-comparing version was built and run: it calls the **correct shipped database
  stale on 116 of 166 markers**. A guard that cries stale on a good database gets switched
  off, so that version is worse than none.
- **The guard's own advice did not work, for the same reason the `--db-root` note above
  exists.** It printed `--pack human --check`, which on a healthy pack exited 3 with "no
  marker database" because of that default. The warning now passes `--markers-dir <the
  directory it judged>`, and that printed command was run verbatim to confirm it. D51
  then removed the underlying trap by deriving the default from `$BTF_DB`, which is the
  difference between a documented workaround and a fixed default.

Verified against the retained pre-D40 database (`backups/markers-pre-D40-rebuild-20261021`,
163 rows): against the current 166-row catalogue it names exactly `salm_invA`, `salm_ttrC`,
`salm_ssaV` and exits 3, and the launch guard run under Nextflow names the same three.

## A behaviour change to be aware of: unevaluable vaccine and veto checks cap the tier

Found by cross-checking all 58 rule sets against the built manifest: 8 of 9 `vaccine_if`
rules and 6 `veto_if` rules referenced markers that can only ever be UNKNOWN. Under Kleene
logic UNKNOWN in those two rule types means *no downgrade* — the unsafe direction — so a
Sabin vaccine strain satisfied Enterovirus C's `confirm_requires` and reached
`confirmed_in_silico`: a vaccine strain reported as confirmed poliovirus.

An unevaluable vaccine or veto check now holds the finding at `supported`. **This reduces
confirmations**, deliberately: "supported, vaccine status untestable" costs a follow-up,
"confirmed wild-type poliovirus" costs an outbreak response. Building the 9 `motif_todo` and
4 `curation_required` markers is what restores those confirmations.

**Re-measured 2026-10-21** against 59 rule sets and the rebuilt manifest, using
`marker_rules.compile_rule` rather than a hand grep: **9 of 10 `vaccine_if`** and **4 of 10
`veto_if`** rules still reference one of the 16 unevaluable markers. `veto_if` is down from 6
because D30 repaired `bthai_ara`; `vaccine_if` went 8 of 9 to 9 of 10 because a rule set was
added, not because one regressed. `confirm_requires` (54 rules) and `support_requires` (31)
are **fully backed** — zero unevaluable references. The four remaining vetoes are
*F. tularensis* (`fnov_specific`), *Variola* (`vaccinia_specific`) and Nipah/Hendra (both
`henipa_cedar`).

## Version pins

All `modules/local/*.nf` conda pins were written from memory and 11 of 13 were stale. Now
updated. Two deliberate exceptions:

- **`kraken2=2.17.2`** — not hygiene. It must match the version that *builds* the index, or
  the format may not read back.
- **`python=3.11` left alone.** Newest is 3.15.0rc2; pinning a pipeline to a release
  candidate is worse than being a version behind.

Several are major-version jumps (fastp 0.23→1.3, chopper 0.7→0.14, sylph 0.6→1.0, mmseqs2
15→18, metabuli 1.0→1.2).

**This warning cashed in.** `samtools=1.24 conda-forge::python=3.11` is **unsatisfiable** —
samtools 1.24 wants `htslib >=1.24,<1.25`, which wants `libdeflate >=1.25,<1.27`, which does
not exist in the channels. It solves happily alongside minimap2, so the conflict is
specifically with a conda-forge python, and MARKER_EVIDENCE died on it during the first
end-to-end run. All 13 distinct specs were then put through the solver: 12 solved, that one
did not. samtools is now pinned to **1.22** everywhere, verified to solve both alongside
minimap2 and alongside python 3.11.

Re-validate with `bash scripts/check_conda_specs.sh` **after any pin change**. It is not in
`tests/test_pipeline.py` because it needs the network and takes minutes; the unit suite stays
offline and fast by design. CI runs it weekly and on demand (`workflow_dispatch`) for the
same reason.

## Environment facts worth not rediscovering

**HOST-SPECIFIC**, all of it — but each of these presented as a bug in something else, which
is why they are written down.

- **rsync (port 873) is blocked** on this host; HTTPS and FTP are open. `kraken2-build`
  defaults to rsync and reports the failure as `"Error downloading assembly summary file"`,
  which looks like an NCBI outage. `--transport auto` probes and falls back.
- **`systemd-oomd` kills anything memory-hungry that runs inside the user session.**
  `/usr/lib/systemd/system/user@.service.d/10-oomd-user-service-defaults.conf` sets
  `ManagedOOMMemoryPressure=kill` with a **50% pressure limit** on `user@1000.service`.
  Bracken loading the 123 GB hash pushed pressure to 87.72%, and oomd killed the build,
  gnome-shell and five tmux scopes at once — on a 503 GB box with no kernel OOM. Pressure,
  not exhaustion, is the trigger. Run long/large jobs as a SYSTEM unit, which is not
  monitored (`ManagedOOMMemoryPressure=auto`):
  `sudo systemd-run --unit=<name> --uid=1000 --gid=1000 /bin/bash <script>`
  (`--uid` keeps the output files owned by the invoking user). Note
  `ManagedOOMMemoryPressure=off` is NOT a valid value — only `auto` or `kill` — so passing
  `off` fails the whole call with `Invalid argument`.
- **A detached background task does not survive the session**, and `setsid` is not enough
  when the launching process is a child subreaper, because an orphan reparents to it rather
  than to PID 1. Same fix: a systemd unit. Verify with `PPID 1`.
- NCBI Entrez is queried **anonymously** (3 req/s). Supplying an email address or
  `NCBI_API_KEY` raises it to 10/s. No address is placed in outbound API URLs by default,
  deliberately: that is the operator's to opt into.
- Nextflow self-updated to **26.04.6**, which enforces strict config/DSL syntax. Any pre-24
  Nextflow pipeline on this machine will now fail to parse.

## Operational work still open

`docs/roadmap.md` is authoritative for priority and for everything not yet built. What
follows is the operational residue — the parts that need a machine, a database or a
measurement rather than a design decision.

1. **Regenerate the ONT plumbing control.** The Illumina positive control exists at
   `$BTF_DB/simulated/plumbing/` — 2M reads, 5000 rpm *Y. pestis*, 94.6% host, truth row in
   `manifest.tsv`. Its ONT partner was cut off by a 1700 s timeout and was never rebuilt.
2. **ONT data on the chemistry the floor was calibrated for.** The only ONT data this
   project has processed is R9.4.1 (see **Real data**), and D10's 0.90 identity floor was
   reasoned for R10.4.1. On R9.4.1 that floor is the only failed gate on every ONT call, so
   the arm's top tier is untested rather than measured — and the floor was not falsified by
   this benchmark, which cannot test it.
3. **Set thresholds from measurement.** Three are explicitly PROVISIONAL and flagged in the
   code: the ONT read-identity floors (D10), `--ont-identity-slack` in the marker tier, and
   the read/breadth gate incoherence (D11) where **75 of 89** taxa carry a decorative
   `min_unique_reads` — up to 243× for *Schistosoma mansoni*. Run
   `scripts/check_threshold_coherence.py` for the table. D12's remedy is also still
   unchosen; the calibration for it is under **Still unresolved**.
4. **A containerised deployment** — roadmap §4, and the project's own standard: required
   before any result is defensible or publishable.
5. **Apply D43's 32-reference cover** for the 9 bacterial species that do not span their own
   diversity. D43 is the measurement and the recommendation; nothing was added. Doing it
   touches nine competitive bundles, which is why D45's rebuild predicate was closed first.

## Real data — first measurement, 2026-10-08

`real/mock_v1`: ZymoBIOMICS D6300 (EVEN + LOG, Illumina and ONT) plus all ten member
organisms as single-organism Illumina libraries. 54 GB, md5-verified. Run as `zymo_rl100`,
`zymo_rl150` and `zymo_ont` because Bracken takes one read length per run.

It refuted D12's neighbour-credibility threshold (now defaulted off) and answered D14's
VERIFY on the tier ladder, cutting real-data false positives 12 → 2 with both simulated
benchmarks unchanged. See **D18**.

**`zymo_ont` COMPLETED 2026-10-16** (15h 36m, 425 CPU-hours) — the first ONT data this
project processed end to end. `MARKER_SEARCH_AA` took **7h 20m** on one sample, so the old
4 h default was off by nearly a factor of two and the 24 h raise was the right size.
Specificity is clean: 84 absent taxa, four findings, all held at `attributed_to_neighbour`,
zero false confirmations. Sensitivity is one of four: *L. monocytogenes* `supported`, the
other three absorbed by `attributed_to_neighbour`.

**Every ONT call is capped by the D10 identity floor and nothing else.** *Listeria* has
breadth 0.95, rpm 136,951, neighbour_ratio 9.7 and BOTH confirmatory markers PRESENT — its
only failed gate is `identity 0.8579 >= 0.9000 (ont/default)`. At 0.90 no ONT sample can
reach `confirmed_in_silico` at all. **Do not lower the floor on this evidence:** the data is
Nicholls 2019, which is R9.4.1, and D10's 0.90 was reasoned for R10.4.1 at ~97-98%. The
floor was not falsified; this benchmark cannot test it. The ONT arm still has no data on the
chemistry it is calibrated for. See **D35**.

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
