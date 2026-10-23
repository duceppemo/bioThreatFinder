#!/usr/bin/env python3
"""measure_reference_span.py - would a known isolate of a target species be
confirmed against the references we ship, or silently degraded?

WHY THIS IS A SCRIPT

D38 measured a genuine *L. monocytogenes* isolate at read identity 0.9643 against the
0.98 Illumina floor and refused to lower the floor. D42 then gave that species a second
reference and the same reads went to 0.9947. Both numbers are real, and the way they were
obtained has a defect that no amount of care in the measurement can fix: THE SPECIES WAS
CHOSEN BECAUSE A SAMPLE COMPLAINED, AND THE IMPROVEMENT WAS SCORED ON THE SAME READS THAT
MOTIVATED IT. No isolate was held out. 30 other target bacteria were never examined,
because no specimen had yet arrived to accuse them.

That is backwards for a surveillance pack. A true sample is an isolate nobody has seen, and
the pack has to be ready for it before it arrives. The question this script answers takes
no specimen at all:

    for every assembly of this species that is publicly known, would reads from THAT
    isolate clear the pipeline's gates against the reference set we ship?

The worst answer over all known assemblies is the pack's guarantee. It is a property of the
pack, not of a run, so it can be measured once, recorded, and re-checked in CI when the
reference set or the thresholds change.

THE SELECTION DEFECT THIS ALSO EXPOSES

build_pack_refs.py ranks candidate assemblies by NCBI metadata alone - refseq_category,
then assembly level, then GCF over GCA, then size (Candidate.sort_key). Nothing in that key
knows what a lineage is, so `--max-per-taxon 2` can return two assemblies from the same
clade and leave the other half of the species uncovered. D42's second reference happened to
land in the other lineage; that was luck in the pick, not a property of the picker. The
`cover` stage here replaces the luck with a stated criterion: choose references so that
every known assembly is within reach of one of them.

THE CALIBRATION TRAP - WHY THE VERDICT IS NOT AN ANI NUMBER

Mash ANI is cheap enough to run over every known assembly, so it is used HERE ONLY TO
CHOOSE WHICH GENOMES TO TEST and to count how many sit near the edge. It is NOT the
verdict, because the conversion from ANI to the pipeline's `mean_identity` is
species-dependent, and not by a little:

    L. monocytogenes   mean_identity crosses 0.98 at mash-ANI ~0.974
    L. pneumophila     mean_identity crosses 0.98 at mash-ANI ~0.981

Read identity is measured over the regions that still align, so how far it floats above
genome-wide ANI depends on how a species' divergence is distributed - uniform drift versus
a few replaced loci. A single ANI cut-off applied across the pack would therefore
under-call some species and over-call others. So the `verdict` stage takes the assemblies
ANI nominated as worst, simulates reads from them, and pushes those reads through THE
PIPELINE'S OWN aligner, THE PIPELINE'S OWN competitive bundle and THE PIPELINE'S OWN
coverage_stats.py. The number that decides is the number the pipeline would itself produce.

WHAT A FAILURE LOOKS LIKE, AND WHY IT IS THE PROJECT'S SIGNATURE BUG

A divergent isolate does not announce itself. Its reads align at reduced identity over
reduced breadth, the identity gate declines to confirm, and the finding settles at
`supported` or lower - which reads exactly like a cautious call on a weak signal rather
than like a reference gap. Worse, at the far end the reads do not align at all: this script
found 9 *C. botulinum* assemblies sharing under 100 of 10,000 hashes with the shipped
Group I reference, i.e. no measurable similarity. Reads from those isolates produce
near-zero breadth, which is indistinguishable from a clean negative. That is the bug shape
this project exists to hunt, sitting in the reference set rather than in the code.

STAGES

    frame    ask NCBI which assemblies of each target species exist (metadata only)
    fetch    download them (bulk; tens of GB for the whole pack)
    span     sketch every assembly, measure its ANI to the shipped reference set,
             and report the minimum, the quantiles and how many sit near the edge
    cover    greedy set cover, SEEDED WITH THE REFERENCES ALREADY SHIPPED, giving the
             smallest set of additions that brings every known assembly within reach
    verdict  simulate reads from the worst assemblies and score them with the real
             pipeline, which is what actually decides

    python3 scripts/measure_reference_span.py --stage frame,fetch,span,cover \\
        --pack-refs $BTF_DB/packs/human/refs \\
        --work $BTF_DB/cache/ani_frame --group bacteria

    python3 scripts/measure_reference_span.py --stage verdict \\
        --pack-refs ... --work ... --out docs/reference_span.tsv

EVERY EXTERNAL TOOL IS RUN BEFORE IT IS USED. D41 records a rebuild that degraded a third
of the marker catalogue because `hmmbuild` resolved to a path on disk and then exited 127
on a missing shared library; the dry run printed the path without ever invoking it.
`require_tool()` here executes each binary and refuses to start if it does not come back.
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import random
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

# distance 0.026 == ANI 0.974; the lower of the two calibrated crossings, used only
# to build the covering graph and to flag edge cases, never as a verdict.
DEFAULT_ANI_FLOOR = 0.974
EDGE_BAND = 0.004
MASH_K = 21
MASH_S = 10000
SIM_COVERAGE = 15
SIM_READ_LEN = 150


def die(msg: str) -> None:
    sys.exit(f"ERROR: {msg}")


def say(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def require_tool(path: str, probe: list[str], label: str) -> str:
    """Run the tool. A resolvable path is not a working tool (D41)."""
    exe = shutil.which(path) or path
    if not Path(exe).exists():
        die(f"{label}: {path!r} not found")
    try:
        proc = subprocess.run([exe] + probe, capture_output=True, text=True, timeout=120)
    except OSError as exc:
        die(f"{label}: {exe} would not execute: {exc}")
    except subprocess.TimeoutExpired:
        die(f"{label}: {exe} hung on {' '.join(probe)}")
    # art_illumina and mash exit non-zero on a bare --version/-h; a tool that
    # produced output ran. Only silence plus failure means broken.
    if proc.returncode != 0 and not (proc.stdout or proc.stderr):
        die(f"{label}: {exe} exited {proc.returncode} with no output; it does not run")
    if "error while loading shared libraries" in (proc.stderr or ""):
        die(f"{label}: {exe} is missing a shared library: {proc.stderr.strip()[:200]}")
    return exe


# ---------------------------------------------------------------------------
# the pack: which taxids are targets and what we ship for them
# ---------------------------------------------------------------------------

def load_pack(refs_dir: Path, group: str | None) -> tuple[dict, dict, dict]:
    """taxid -> [shipped accessions], taxid -> name, taxid -> group."""
    man = refs_dir / "manifest.tsv"
    if not man.exists():
        die(f"no manifest at {man}")
    shipped: dict = collections.defaultdict(list)
    names: dict = {}
    groups: dict = {}
    with man.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_tax, i_name, i_grp = (header.index("taxid"), header.index("name"),
                                    header.index("group"))
            i_role, i_acc = header.index("role"), header.index("accession")
        except ValueError:
            die(f"{man}: unexpected header {header}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_tax, i_name, i_grp, i_role, i_acc):
                continue
            if f[i_role] != "target":
                continue
            if group and f[i_grp] != group:
                continue
            if not f[i_acc]:
                continue
            shipped[f[i_tax]].append(f[i_acc])
            names[f[i_tax]] = f[i_name]
            groups[f[i_tax]] = f[i_grp]
    if not shipped:
        die(f"{man}: no target rows" + (f" in group {group}" if group else ""))
    return shipped, names, groups


# ---------------------------------------------------------------------------
# frame: what exists
# ---------------------------------------------------------------------------

def frame_path(work: Path, taxid: str) -> Path:
    return work / taxid / "frame.tsv"


def stage_frame(taxids, names, work: Path, datasets: str, dataformat: str,
                levels: str, force: bool) -> None:
    """The sampling frame is the taxid SUBTREE at RefSeq complete/chromosome level.

    A subtree is not a species. Several target taxids carry subspecies that are
    arguably separate organisms - *C. botulinum* spans four phenotypically-defined
    groups with no close genomic relationship, and taxid 263 brings in novicida-like
    entries - so a very low ANI inside a frame can be real biology, a taxonomy
    artefact, or both. The frame keeps the metadata so the report can say which.
    """
    for taxid in taxids:
        out = frame_path(work, taxid)
        if out.exists() and out.stat().st_size and not force:
            say(f"  frame {taxid} {names.get(taxid,'')}: cached "
                f"({sum(1 for _ in out.open())} assemblies)")
            continue
        out.parent.mkdir(parents=True, exist_ok=True)
        p1 = subprocess.run(
            [datasets, "summary", "genome", "taxon", taxid,
             "--assembly-source", "RefSeq", "--assembly-level", levels,
             "--exclude-atypical", "--mag", "exclude", "--as-json-lines"],
            capture_output=True, text=True)
        if p1.returncode != 0:
            say(f"  frame {taxid}: datasets failed: {p1.stderr.strip()[:200]}")
            continue
        p2 = subprocess.run(
            [dataformat, "tsv", "genome", "--fields",
             "accession,assminfo-level,organism-name,"
             "organism-infraspecific-strain,organism-infraspecific-isolate,"
             "assmstats-total-sequence-len,assmstats-number-of-contigs"],
            input=p1.stdout, capture_output=True, text=True)
        if p2.returncode != 0:
            say(f"  frame {taxid}: dataformat failed: {p2.stderr.strip()[:200]}")
            continue
        rows = [ln for ln in p2.stdout.splitlines()[1:] if ln.strip()]
        out.write_text("\n".join(rows) + ("\n" if rows else ""))
        say(f"  frame {taxid} {names.get(taxid,'')}: {len(rows)} assemblies")


def read_frame(work: Path, taxid: str) -> dict:
    """accession -> (level, organism, strain, isolate, length, n_contigs)."""
    p = frame_path(work, taxid)
    if not p.exists():
        return {}
    out = {}
    for line in p.open():
        f = line.rstrip("\n").split("\t")
        if not f or not f[0]:
            continue
        def g(i):
            return f[i] if len(f) > i else ""
        try:
            ncont = int(g(6))
        except ValueError:
            ncont = 99
        out[f[0]] = (g(1), g(2), g(3), g(4), g(5), ncont)
    return out


# ---------------------------------------------------------------------------
# fetch: bulk download, into its own tree
# ---------------------------------------------------------------------------

def stage_fetch(taxids, names, work: Path, datasets: str, batch: int) -> None:
    for taxid in taxids:
        frame = read_frame(work, taxid)
        gen = work / taxid / "genomes"
        gen.mkdir(parents=True, exist_ok=True)
        todo = [a for a in sorted(frame) if not (gen / f"{a}.fna").exists()]
        if not todo:
            say(f"  fetch {taxid} {names.get(taxid,'')}: {len(frame)} on disk")
            continue
        say(f"  fetch {taxid} {names.get(taxid,'')}: {len(todo)} of {len(frame)} to get")
        tmp = work / taxid / "tmp"
        tmp.mkdir(parents=True, exist_ok=True)
        for i in range(0, len(todo), batch):
            chunk = todo[i:i + batch]
            acc_file = tmp / "acc.txt"
            acc_file.write_text("\n".join(chunk) + "\n")
            zp = tmp / "p.zip"
            un = tmp / "un"
            shutil.rmtree(un, ignore_errors=True)
            un.mkdir(parents=True, exist_ok=True)
            proc = subprocess.run(
                [datasets, "download", "genome", "accession", "--inputfile",
                 str(acc_file), "--include", "genome", "--no-progressbar",
                 "--filename", str(zp)], capture_output=True, text=True)
            if proc.returncode != 0:
                say(f"    WARN batch at {i} failed: {proc.stderr.strip()[:160]}")
                continue
            subprocess.run(["unzip", "-q", "-o", str(zp), "-d", str(un)],
                           capture_output=True, text=True)
            for fna in un.rglob("*_genomic.fna"):
                shutil.move(str(fna), str(gen / f"{fna.parent.name}.fna"))
            shutil.rmtree(un, ignore_errors=True)
            zp.unlink(missing_ok=True)
            say(f"    {sum(1 for _ in gen.glob('*.fna'))} / {len(frame)} on disk")


# ---------------------------------------------------------------------------
# span: ANI of every known assembly to the shipped reference set
# ---------------------------------------------------------------------------

def stage_span(taxids, shipped, names, work: Path, refs_dir: Path, mash: str,
               force: bool) -> None:
    """Sketch every assembly, then measure each against the SHIPPED FASTA.

    The shipped file is the masked one, because that is what reads are aligned to. A
    masked interval contributes no k-mers, which is correct: it cannot recruit a read
    either.
    """
    for taxid in taxids:
        frame = read_frame(work, taxid)
        gen = work / taxid / "genomes"
        have = sorted(gen.glob("*.fna"))
        if not have:
            say(f"  span {taxid}: nothing fetched")
            continue
        if len(have) < len(frame):
            say(f"  span {taxid} {names.get(taxid,'')}: INCOMPLETE "
                f"({len(have)}/{len(frame)}) - the minimum is provisional")
        mdir = work / taxid / "mash"
        mdir.mkdir(parents=True, exist_ok=True)
        qmsh = mdir / "queries.msh"
        if not qmsh.exists() or force:
            lst = mdir / "qlist.txt"
            lst.write_text("\n".join(str(p) for p in have) + "\n")
            run_mash(mash, ["sketch", "-k", str(MASH_K), "-s", str(MASH_S),
                            "-o", str(mdir / "queries"), "-l", str(lst)])
        out = mdir / "dist_long.tsv"
        if out.exists() and not force:
            say(f"  span {taxid} {names.get(taxid,'')}: cached")
            continue
        rows = []
        for n, acc in enumerate(shipped[taxid]):
            fasta = find_shipped_fasta(refs_dir, taxid, acc, names.get(taxid, ""),
                                       is_first=(n == 0))
            if not fasta:
                say(f"  span {taxid}: shipped {acc} has no FASTA in {refs_dir}")
                continue
            rmsh = mdir / f"{acc}.msh"
            if not rmsh.exists() or force:
                run_mash(mash, ["sketch", "-k", str(MASH_K), "-s", str(MASH_S),
                                "-o", str(mdir / acc), str(fasta)])
            proc = run_mash(mash, ["dist", str(rmsh), str(qmsh)])
            for line in proc.splitlines():
                f = line.split("\t")
                if len(f) < 5:
                    continue
                q = Path(f[1]).name
                if q.endswith(".fna"):
                    q = q[:-4]
                rows.append(f"{q}\t{acc}\t{1.0 - float(f[2]):.6f}\t{f[4]}")
        out.write_text("\n".join(rows) + ("\n" if rows else ""))
        say(f"  span {taxid} {names.get(taxid,'')}: {len(rows)} pairs")


def find_shipped_fasta(refs_dir: Path, taxid: str, acc: str, name: str,
                       is_first: bool) -> Path | None:
    """The pack writes the FIRST reference of a taxid under the plain name and every
    later one with the accession infixed (build_pack_refs.py).

    So the plain file may only be claimed by the first manifest row for that taxid.
    Handing it to a later accession because its own file happened to be missing would
    measure the wrong genome and report the answer as if it were that accession's -
    silently, and in the direction that flatters the pack. Return None instead and let
    the caller complain.
    """
    slug = name.replace(" ", "_")
    infixed = refs_dir / f"{taxid}_{slug}.{acc}.fna.gz"
    if infixed.exists():
        return infixed
    hits = sorted(refs_dir.glob(f"{taxid}_*.{acc}.fna.gz"))
    if hits:
        return hits[0]
    if is_first:
        plain = refs_dir / f"{taxid}_{slug}.fna.gz"
        if plain.exists():
            return plain
        bare = [p for p in sorted(refs_dir.glob(f"{taxid}_*.fna.gz"))
                if p.name.count(".") == 2]      # <taxid>_<slug>.fna.gz only
        if len(bare) == 1:
            return bare[0]
    return None


def run_mash(mash: str, args: list[str]) -> str:
    proc = subprocess.run([mash] + args, capture_output=True, text=True)
    if proc.returncode != 0:
        die(f"mash {args[0]} failed: {proc.stderr.strip()[:300]}")
    return proc.stdout


def best_ani(work: Path, taxid: str) -> tuple[dict, dict]:
    """accession -> best ANI over the shipped set; accession -> shared-hash string."""
    best: dict = {}
    shared: dict = {}
    p = work / taxid / "mash" / "dist_long.tsv"
    if not p.exists():
        return best, shared
    for line in p.open():
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        q, _r, ani, sh = f[0], f[1], float(f[2]), f[3]
        if q not in best or ani > best[q]:
            best[q] = ani
            shared[q] = sh
    return best, shared


# ---------------------------------------------------------------------------
# cover: the smallest set of additions that reaches everything
# ---------------------------------------------------------------------------

def stage_cover(taxids, shipped, names, work: Path, mash: str, floor: float,
                force: bool) -> dict:
    """Greedy set cover, seeded with what the pack already ships.

    Seeding matters: the output is then a list of ADDITIONS, so acting on it never
    churns a reference that other recorded measurements were made against.
    """
    maxdist = f"{1.0 - floor:.6f}"
    result: dict = {}
    for taxid in taxids:
        mdir = work / taxid / "mash"
        qmsh = mdir / "queries.msh"
        if not qmsh.exists():
            continue
        cov_file = mdir / "cover.tsv"
        if not cov_file.exists() or force:
            proc = run_mash(mash, ["dist", "-d", maxdist, str(qmsh), str(qmsh)])
            pairs = []
            for line in proc.splitlines():
                f = line.split("\t")
                if len(f) < 3:
                    continue
                a, b = Path(f[0]).name, Path(f[1]).name
                a = a[:-4] if a.endswith(".fna") else a
                b = b[:-4] if b.endswith(".fna") else b
                pairs.append(f"{a}\t{b}")
            cov_file.write_text("\n".join(pairs) + ("\n" if pairs else ""))
        frame = read_frame(work, taxid)
        universe = set(frame)
        covers: dict = collections.defaultdict(set)
        for line in cov_file.open():
            f = line.rstrip("\n").split("\t")
            if len(f) == 2 and f[0] in universe and f[1] in universe:
                covers[f[0]].add(f[1])
                covers[f[1]].add(f[0])
        for a in universe:
            covers[a].add(a)
        best, _sh = best_ani(work, taxid)
        covered = {a for a in universe if best.get(a, 0.0) >= floor}
        # Ties go to the assembly with fewest contigs: a single-sequence complete
        # genome is a better reference, and score_evidence's per-assembly geometry
        # substitution is exact only when every assembly under the taxid is one
        # sequence (D42).
        added = greedy_cover(universe, covers, covered,
                             {a: frame[a][5] for a in universe})
        reached = set(covered)
        for a, _g in added:
            reached |= covers[a]
        result[taxid] = {"added": added, "uncovered": len(universe) - len(reached),
                         "n": len(universe)}
        say(f"  cover {taxid} {names.get(taxid,'')}: ships {len(shipped[taxid])}, "
            f"needs {len(added)} more"
            + (f", {result[taxid]['uncovered']} uncoverable" if result[taxid]['uncovered'] else ""))
    return result


# ---------------------------------------------------------------------------
# verdict: the pipeline's own answer on the worst assemblies
# ---------------------------------------------------------------------------

def stage_verdict(taxids, names, groups, work: Path, refs_dir: Path, repo: Path,
                  art: str, minimap2: str, samtools: str, python: str,
                  n_worst: int, force: bool) -> dict:
    """Simulate, align with the pipeline's command, score with the pipeline's scorer.

    The read set is the whole simulated library rather than the subset a stage-1
    classifier would have assigned. That changes the breadth denominator not at all and
    identity not at all - identity is computed over aligned reads only - and D42 checked
    it directly: the one-reference figure measured this way reproduced the production
    run's 0.9643 to four places.
    """
    out: dict = {}
    for taxid in taxids:
        best, shared = best_ani(work, taxid)
        if not best:
            continue
        mmi = refs_dir / "competitive" / f"{taxid}.sr.mmi"
        if not mmi.exists():
            say(f"  verdict {taxid}: no competitive bundle, skipped")
            continue
        order = sorted(best.items(), key=lambda kv: kv[1])
        panel = order[:n_worst] + [order[len(order) // 2]]
        seen = set()
        rows = []
        cal = work / taxid / "calib"
        cal.mkdir(parents=True, exist_ok=True)
        for acc, ani in panel:
            if acc in seen:
                continue
            seen.add(acc)
            cov_tsv = cal / f"{acc}.cov.tsv"
            if not cov_tsv.exists() or force:
                simulate_and_score(acc, taxid, groups.get(taxid, "bacteria"),
                                   work, refs_dir, repo, cal, mmi, art, minimap2,
                                   samtools, python)
            rows.append((acc, ani, shared.get(acc, ""), read_cov_row(cov_tsv, taxid)))
        out[taxid] = rows
        for acc, ani, sh, m in rows:
            say(f"  verdict {taxid} {names.get(taxid,'')[:24]:<24} {acc} "
                f"ani={ani:.4f} shared={sh} identity={m.get('mean_identity','?')} "
                f"breadth={m.get('breadth','?')}")
    return out


def simulate_and_score(acc, taxid, group, work: Path, refs_dir: Path, repo: Path,
                       cal: Path, mmi: Path, art, minimap2, samtools, python) -> None:
    genome = work / taxid / "genomes" / f"{acc}.fna"
    if not genome.exists():
        return
    stem = cal / acc
    subprocess.run([art, "-ss", "HS25", "-i", str(genome), "-p", "-l",
                    str(SIM_READ_LEN), "-f", str(SIM_COVERAGE), "-m", "350",
                    "-s", "50", "-na", "-rs", "42", "-o", f"{stem}."],
                   capture_output=True, text=True)
    r1, r2 = Path(f"{stem}.1.fq"), Path(f"{stem}.2.fq")
    if not r1.exists():
        say(f"    {acc}: simulation produced no reads")
        return
    bam = Path(f"{stem}.bam")
    mm = subprocess.Popen([minimap2, "-ax", "sr", "-t", "8", "--secondary=no",
                           str(mmi), str(r1), str(r2)],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    st = subprocess.run([samtools, "sort", "-@", "4", "-o", str(bam), "-"],
                        stdin=mm.stdout, capture_output=True, text=True)
    mm.wait()
    if st.returncode != 0:
        say(f"    {acc}: sort failed: {st.stderr.strip()[:160]}")
        return
    subprocess.run([samtools, "index", str(bam)], capture_output=True, text=True)
    env = dict(os.environ)
    env["PATH"] = f"{Path(samtools).parent}:{env.get('PATH','')}"
    proc = subprocess.run(
        [python, "-I", str(repo / "bin" / "coverage_stats.py"), "--bam", str(bam),
         "--sample", acc, "--taxid", taxid, "--group", group,
         "--platform", "illumina", "--mask-bed", str(refs_dir / "masked"),
         "--out", str(cal / f"{acc}.cov.tsv")],
        capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        say(f"    {acc}: coverage_stats failed: {proc.stderr.strip()[:200]}")
    for junk in (r1, r2, bam, Path(f"{bam}.bai")):
        junk.unlink(missing_ok=True)


def read_cov_row(path: Path, taxid: str) -> dict:
    if not path.exists():
        return {}
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            row = dict(zip(header, f))
            if row.get("taxid") == taxid:
                return row
    return {}


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

def quantile(vals, p):
    v = sorted(vals)
    if not v:
        return float("nan")
    i = max(0, min(len(v) - 1, int(round(p / 100.0 * (len(v) - 1)))))
    return v[i]


def write_report(taxids, shipped, names, groups, work: Path, cover: dict,
                 verdict: dict, floor: float, out: Path | None) -> None:
    cols = ["taxid", "name", "group", "n_known", "n_shipped_refs", "min_ani",
            "p1_ani", "p5_ani", "median_ani", "n_below_floor", "pct_below_floor",
            "n_edge", "n_unmeasurable", "refs_to_add", "add_accessions",
            "worst_accession", "worst_ani", "worst_mean_identity", "worst_breadth",
            "worst_best_ref_breadth", "frame_complete"]
    lines = ["\t".join(cols)]
    for taxid in taxids:
        best, shared = best_ani(work, taxid)
        if not best:
            continue
        frame = read_frame(work, taxid)
        vals = list(best.values())
        nbad = sum(1 for v in vals if v < floor)
        nedge = sum(1 for v in vals if abs(v - floor) <= EDGE_BAND)
        nunmeas = sum(1 for s in shared.values()
                      if s and s.split("/")[0].isdigit() and int(s.split("/")[0]) < 100)
        cv = cover.get(taxid, {})
        added = cv.get("added", [])
        vd = verdict.get(taxid, [])
        w = vd[0] if vd else None
        m = w[3] if w else {}
        lines.append("\t".join(str(x) for x in [
            taxid, names.get(taxid, ""), groups.get(taxid, ""), len(best),
            len(shipped.get(taxid, [])), f"{min(vals):.4f}",
            f"{quantile(vals,1):.4f}", f"{quantile(vals,5):.4f}",
            f"{quantile(vals,50):.4f}", nbad, f"{100.0*nbad/len(vals):.1f}",
            nedge, nunmeas, len(added),
            ";".join(a for a, _g in added) or "-",
            w[0] if w else "-", f"{w[1]:.4f}" if w else "-",
            m.get("mean_identity", "-"), m.get("breadth", "-"),
            m.get("best_reference_breadth", "-"),
            "yes" if len(best) >= len(frame) else "no",
        ]))
    text = "\n".join(lines) + "\n"
    if out:
        out.write_text(text)
        say(f"wrote {out}")
    else:
        print(text, end="")


def greedy_cover(universe, covers, covered, contigs) -> list:
    """Shared by the cover stage and the self-test, so the test exercises the real thing."""
    added = []
    covered = set(covered)
    while len(covered) < len(universe):
        pick = max(universe, key=lambda s: (len(covers[s] - covered),
                                            -contigs.get(s, 99), s))
        gain = len(covers[pick] - covered)
        if gain == 0:
            break
        covered |= covers[pick]
        added.append((pick, gain))
    return added


def self_test() -> int:
    """Offline tests for the three things that decide what gets REPORTED.

    None needs a pack, a genome or a network. Each was, or could be, wrong in a way
    that yields a confident-looking table rather than an error - which for a script
    whose whole purpose is to find silent degradation would be the same bug one level up.
    """
    failures: list[str] = []

    def check(name: str, cond: bool, detail: str = "") -> None:
        print(f"  {'ok  ' if cond else 'FAIL'} {name}"
              + (f": {detail}" if detail and not cond else ""))
        if not cond:
            failures.append(name)

    print("measure_reference_span.py --self-test")
    print("-- the plain filename belongs to the first shipped reference, nobody else --")

    with tempfile.TemporaryDirectory(prefix="refspan.") as td:
        refs = Path(td)
        (refs / "1639_Listeria_monocytogenes.fna.gz").write_bytes(b"")
        (refs / "1639_Listeria_monocytogenes.GCF_000008285.1.fna.gz").write_bytes(b"")
        first = find_shipped_fasta(refs, "1639", "GCF_000196035.1",
                                   "Listeria monocytogenes", is_first=True)
        check("the first accession resolves to the plain name",
              first is not None and first.name == "1639_Listeria_monocytogenes.fna.gz",
              str(first))
        second = find_shipped_fasta(refs, "1639", "GCF_000008285.1",
                                    "Listeria monocytogenes", is_first=False)
        check("a later accession resolves to its own infixed file",
              second is not None and second.name.endswith(".GCF_000008285.1.fna.gz"),
              str(second))
        # The one that matters: a later accession whose file is absent must NOT be
        # handed the plain file. That would measure the first genome and report the
        # number under the second accession's name - wrong, silent, and flattering,
        # because the plain file is the better-curated assembly.
        orphan = find_shipped_fasta(refs, "1639", "GCF_999999999.9",
                                    "Listeria monocytogenes", is_first=False)
        check("a later accession with no file of its own gets None, not the plain file",
              orphan is None, f"returned {orphan}")

    print("-- the cover is seeded, so it reports additions and never churn --")
    # Three clusters of two. Shipped reference covers cluster A only.
    universe = {"a1", "a2", "b1", "b2", "c1", "c2"}
    covers = {"a1": {"a1", "a2"}, "a2": {"a1", "a2"},
              "b1": {"b1", "b2"}, "b2": {"b1", "b2"},
              "c1": {"c1", "c2"}, "c2": {"c1", "c2"}}
    contigs = {"a1": 1, "a2": 1, "b1": 7, "b2": 1, "c1": 1, "c2": 3}
    added = greedy_cover(universe, covers, {"a1", "a2"}, contigs)
    check("two uncovered clusters need exactly two additions",
          len(added) == 2, f"{added}")
    check("and neither addition comes from the cluster already covered",
          all(a[0] not in ("a1", "a2") for a in added), f"{added}")
    check("a tie inside a cluster goes to the single-sequence assembly",
          "b2" in [a for a, _g in added] and "b1" not in [a for a, _g in added],
          f"{added}")
    already = greedy_cover(universe, covers, set(universe), contigs)
    check("a pack that already covers everything needs no additions",
          already == [], f"{already}")

    print("-- a near-zero shared-hash count is not a low ANI, it is no measurement --")
    # mash reports a distance for 11/10000 shared hashes, and 1-d looks like an ANI of
    # 0.7086. It is noise: below ~100 shared hashes there is nothing to estimate from.
    # Counting those separately is what stops the report claiming a precise ANI for a
    # genome its reference cannot see at all.
    shared = {"g1": "11/10000", "g2": "4780/10000", "g3": "99/10000",
              "g4": "100/10000"}
    nunmeas = sum(1 for s in shared.values()
                  if s and s.split("/")[0].isdigit() and int(s.split("/")[0]) < 100)
    check("shared-hash counts below 100 are flagged unmeasurable", nunmeas == 2,
          f"got {nunmeas}")

    print("-- a tool that resolves but does not run is not a usable tool (D41) --")
    with tempfile.TemporaryDirectory(prefix="refspan.") as td:
        broken = Path(td) / "faketool"
        broken.write_text("#!/bin/sh\necho 'faketool: error while loading shared "
                          "libraries: libnope.so.1' >&2\nexit 127\n")
        broken.chmod(0o755)
        rc = subprocess.run([sys.executable, __file__, "--self-test-probe",
                             str(broken)], capture_output=True, text=True)
        check("require_tool refuses a binary that cannot load its libraries",
              rc.returncode != 0 and "shared library" in (rc.stdout + rc.stderr),
              f"rc={rc.returncode} out={(rc.stdout + rc.stderr).strip()[:120]}")
        good = Path(td) / "goodtool"
        good.write_text("#!/bin/sh\necho ok\n")
        good.chmod(0o755)
        rc = subprocess.run([sys.executable, __file__, "--self-test-probe",
                             str(good)], capture_output=True, text=True)
        check("and accepts one that runs", rc.returncode == 0,
              f"rc={rc.returncode} out={(rc.stdout + rc.stderr).strip()[:120]}")

    if failures:
        print(f"\n{len(failures)} self-test(s) FAILED: {', '.join(failures)}")
        return 1
    print("all self-tests passed")
    return 0


def main() -> int:
    # Tiny helper process used by the self-test to prove require_tool rejects a
    # binary that resolves and then fails to execute.
    if len(sys.argv) == 3 and sys.argv[1] == "--self-test-probe":
        require_tool(sys.argv[2], ["--version"], "probe")
        return 0

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--self-test", action="store_true",
                    help="run offline unit tests and exit; needs no pack and no network")
    ap.add_argument("--pack-refs", required=False, type=Path,
                    help="the pack's refs directory, holding manifest.tsv, "
                         "competitive/ and masked/")
    ap.add_argument("--work", required=False, type=Path,
                    help="where frames, downloaded assemblies and sketches live. "
                         "Tens of GB for the whole bacterial pack; put it on the "
                         "data volume, not in the repo.")
    ap.add_argument("--stage", default="span,cover",
                    help="comma-separated: frame,fetch,span,cover,verdict")
    ap.add_argument("--group", default="bacteria",
                    help="restrict to one manifest group; empty means all. The "
                         "default is bacteria because the 0.98 identity floor this "
                         "measures against is the bacterial one - a divergent virus "
                         "is legitimately far from any reference and has its own "
                         "floor (D10).")
    ap.add_argument("--taxid", action="append", default=[],
                    help="restrict to these taxids; repeatable")
    ap.add_argument("--ani-floor", type=float, default=DEFAULT_ANI_FLOOR,
                    help=f"ANI used to build the covering graph and flag edge cases "
                         f"[{DEFAULT_ANI_FLOOR}]. NOT a verdict: the ANI at which "
                         f"mean_identity crosses 0.98 is species-dependent (0.974 for "
                         f"L. monocytogenes, 0.981 for L. pneumophila), which is why "
                         f"the verdict stage simulates instead.")
    ap.add_argument("--assembly-level", default="complete,chromosome",
                    help="frame stage: which assembly levels count as known")
    ap.add_argument("--n-worst", type=int, default=2,
                    help="verdict stage: how many of the furthest assemblies to "
                         "simulate per species, plus the median one [2]")
    ap.add_argument("--batch", type=int, default=200,
                    help="fetch stage: accessions per download request")
    ap.add_argument("--out", type=Path, default=None, help="report TSV")
    ap.add_argument("--force", action="store_true",
                    help="recompute cached stages. Existence is not freshness (D42).")
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    ap.add_argument("--datasets", default="datasets")
    ap.add_argument("--dataformat", default="dataformat")
    ap.add_argument("--mash", default="mash")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--art", default="art_illumina")
    ap.add_argument("--python", default=sys.executable)
    args = ap.parse_args()

    if args.self_test:
        return self_test()
    # NOT required=True above: --self-test needs neither a pack nor a work tree, and
    # argparse enforces required arguments before it can be seen.
    missing = [n for n, v in (("--pack-refs", args.pack_refs), ("--work", args.work))
               if v is None]
    if missing:
        ap.error(f"{' and '.join(missing)} required (except with --self-test)")

    stages = [s.strip() for s in args.stage.split(",") if s.strip()]
    bad = [s for s in stages if s not in ("frame", "fetch", "span", "cover", "verdict")]
    if bad:
        die(f"unknown stage(s): {', '.join(bad)}")

    shipped, names, groups = load_pack(args.pack_refs, args.group or None)
    taxids = sorted(shipped, key=lambda t: int(t) if t.isdigit() else 0)
    if args.taxid:
        want = set(args.taxid)
        taxids = [t for t in taxids if t in want]
        if not taxids:
            die("no target taxid matched --taxid")
    args.work.mkdir(parents=True, exist_ok=True)
    say(f"{len(taxids)} target taxid(s); stages: {', '.join(stages)}")

    if "frame" in stages:
        ds = require_tool(args.datasets, ["--version"], "datasets")
        df = require_tool(args.dataformat, ["--help"], "dataformat")
        stage_frame(taxids, names, args.work, ds, df, args.assembly_level, args.force)
    if "fetch" in stages:
        ds = require_tool(args.datasets, ["--version"], "datasets")
        require_tool("unzip", ["-v"], "unzip")
        stage_fetch(taxids, names, args.work, ds, args.batch)
    cover: dict = {}
    verdict: dict = {}
    if "span" in stages:
        mash = require_tool(args.mash, ["--version"], "mash")
        stage_span(taxids, shipped, names, args.work, args.pack_refs, mash, args.force)
    if "cover" in stages:
        mash = require_tool(args.mash, ["--version"], "mash")
        cover = stage_cover(taxids, shipped, names, args.work, mash,
                            args.ani_floor, args.force)
    if "verdict" in stages:
        art = require_tool(args.art, ["-h"], "art_illumina")
        mm2 = require_tool(args.minimap2, ["--version"], "minimap2")
        st = require_tool(args.samtools, ["--version"], "samtools")
        py = require_tool(args.python, ["-V"], "python")
        verdict = stage_verdict(taxids, names, groups, args.work, args.pack_refs,
                                args.repo, art, mm2, st, py, args.n_worst, args.force)
    write_report(taxids, shipped, names, groups, args.work, cover, verdict,
                 args.ani_floor, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
