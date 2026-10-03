#!/usr/bin/env python3
"""fetch_mock_communities.py - download the real mock-community runs and emit a
benchmark manifest the evaluator can read.

WHY THIS EXISTS AS A SCRIPT RATHER THAN A FEW curl LINES

Every failure mode in this fetch is SILENT. A truncated download yields fewer reads,
which reads downstream as "fewer organisms detected"; a wrong accession measures a
different library and reports it under this one's name; the wrong Bracken read-length
model produces plausible abundances that are simply incorrect. None of these error.
So each one gets an explicit check that fails loudly instead.

Specifically:

  * md5 is verified against ENA for every file. Byte count alone is not enough - a
    resumed partial download can match a stale expected size.
  * runs with read_count == 0 are REFUSED. ENA study PRJEB29504 contains six such runs
    (ERR2887847-52) carrying correct-looking titles and platforms. Picking one by title
    gives an empty FASTQ that looks exactly like a clean negative.
  * the ACTUAL read length is measured from the first FASTQ record and compared with the
    scenario's declared `bracken_readlen`. base_count/read_count from ENA is only a hint:
    read_count means spots in some submissions and individual reads in others, so the
    same ratio can mean 2x75 or 1x150.
  * accession metadata is re-queried at fetch time rather than hardcoded, so the URL,
    md5 and size cannot go stale relative to the archive.

The manifest it writes carries the columns bin/evaluate_benchmark.py requires -
sample, scenario_id, platform, truth_present, truth_absent - so real-data results score
through exactly the same path as the simulated batches. Truth comes from the scenario
file, because for a defined community the manufacturer's composition IS the truth; there
is no spike plan to read it from.

NOT A SUBSTITUTE FOR THE SIMULATED SERIES. These standards contain no select agents, so
LoD for the high-consequence watchlist still comes from benchmarks/scenarios/human_core.yaml.

Usage:
  fetch_mock_communities.py --scenarios benchmarks/scenarios/real_mock_communities.yaml \\
                            --out-dir /media/8tb/bioThreatFinder/real/zymo_v1
  fetch_mock_communities.py ... --only zymo_even --subsample 2000000
  fetch_mock_communities.py ... --dry-run
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS = ("run_accession,instrument_platform,instrument_model,library_strategy,"
              "scientific_name,sample_title,tax_id,read_count,base_count,"
              "fastq_ftp,fastq_md5,fastq_bytes")

# Runs that exist, look right and contain nothing. Hardcoded because a guard that only
# checks read_count would still let a future empty run through with a non-zero count.
KNOWN_EMPTY = {"ERR2887847", "ERR2887848", "ERR2887849",
               "ERR2887850", "ERR2887851", "ERR2887852"}


def log(msg: str) -> None:
    print(f"[fetch] {msg}", file=sys.stderr, flush=True)


def die(msg: str) -> None:
    print(f"[fetch] ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


def ena_runs(accession: str) -> list[dict]:
    """Query ENA for a run or study. Raises rather than returning junk."""
    url = f"{ENA_API}?accession={accession}&result=read_run&fields={ENA_FIELDS}&format=tsv"
    try:
        with urllib.request.urlopen(url, timeout=120) as fh:
            text = fh.read().decode()
    except Exception as exc:                                   # noqa: BLE001
        die(f"ENA query failed for {accession}: {exc}")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if not rows:
        die(f"ENA returned no run for {accession}. A typo here would otherwise be "
            f"indistinguishable from an empty download.")
    return rows


def md5_of(path: Path) -> str:
    h = hashlib.md5()                                          # noqa: S324 - ENA publishes md5
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def first_read_length(fq: Path) -> int:
    """Length of the first sequence record, measured rather than inferred."""
    opener = gzip.open if fq.suffix == ".gz" else open
    with opener(fq, "rt") as fh:                               # type: ignore[operator]
        fh.readline()
        return len(fh.readline().strip())


def free_gib(path: Path) -> float:
    probe = path
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    return shutil.disk_usage(probe).free / 2 ** 30


def download(url: str, dest: Path, want_md5: str, dry: bool) -> None:
    stamp = dest.with_suffix(dest.suffix + ".md5ok")
    if stamp.exists() and dest.exists():
        log(f"SKIP {dest.name} (md5 already verified)")
        return
    if dry:
        log(f"[dry-run] would fetch {url} -> {dest}")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    log(f"fetching {dest.name}")
    # curl over https: rsync port 873 is firewalled on this host, see docs/STATUS.md.
    # -C - resumes, which is why md5 rather than size is the acceptance test.
    rc = subprocess.call(["curl", "-sSL", "--fail", "-C", "-", "-o", str(dest),
                          f"https://{url}" if not url.startswith("http") else url])
    if rc != 0:
        die(f"curl exited {rc} for {url}. Partial file left at {dest} for inspection.")
    got = md5_of(dest)
    if got != want_md5:
        die(f"md5 mismatch for {dest.name}: got {got}, ENA says {want_md5}. A truncated "
            f"or corrupted FASTQ reads downstream as fewer organisms, not as an error. "
            f"Delete the file and retry.")
    stamp.write_text(f"{got}\n")
    log(f"  md5 OK {dest.name}")


def fetch_accession(acc: str, out_dir: Path, declared_len: int | None,
                    dry: bool) -> dict:
    if acc in KNOWN_EMPTY:
        die(f"{acc} is one of the six EMPTY runs in PRJEB29504 (read_count=0, no files). "
            f"It carries a correct-looking title and platform, so using it would look "
            f"like a clean negative. Pick a different run.")
    row = ena_runs(acc)[0]
    n_reads = int(row.get("read_count") or 0)
    if n_reads == 0:
        die(f"{acc} reports read_count=0 at ENA, so there is nothing to download. "
            f"Refusing rather than producing an empty FASTQ.")

    ftps = [u for u in (row.get("fastq_ftp") or "").split(";") if u]
    md5s = [m for m in (row.get("fastq_md5") or "").split(";") if m]
    if not ftps:
        die(f"{acc} has no fastq_ftp at ENA.")
    if len(ftps) != len(md5s):
        die(f"{acc}: {len(ftps)} fastq url(s) but {len(md5s)} md5(s) - cannot verify.")

    log(f"{acc}: {row.get('sample_title')} | {row.get('instrument_platform')} "
        f"{row.get('instrument_model')} | tax_id={row.get('tax_id')} | "
        f"{n_reads:,} reads | {len(ftps)} file(s)")

    local = []
    for url, m in zip(ftps, md5s):
        dest = out_dir / "reads" / Path(url).name
        download(url, dest, m, dry)
        local.append(dest)

    if not dry and declared_len is not None:
        actual = first_read_length(local[0])
        hint = (int(row["base_count"]) / n_reads) if n_reads else 0
        log(f"  read length: measured {actual} bp (ENA ratio hint {hint:.0f} bases/spot)")
        if abs(actual - declared_len) > 2:
            die(f"{acc}: measured read length {actual} bp does not match the scenario's "
                f"declared bracken_readlen {declared_len}. Bracken's k-mer distribution "
                f"is specific to the read length, so proceeding would produce plausible "
                f"but WRONG abundances - see nextflow.config:174. Fix bracken_readlen in "
                f"the scenario, or choose the matching distribution.")
    return {"accession": acc, "row": row, "files": local}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenarios", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--only", action="append", default=[],
                    help="scenario id; repeatable. Default: all with accessions.")
    ap.add_argument("--min-free-gib", type=float, default=80.0)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))
    import btf_yaml                                            # noqa: PLC0415

    doc = btf_yaml.load(args.scenarios)
    scenarios = {s["id"]: s for s in (doc.get("scenarios") or [])}
    if not scenarios:
        die(f"no scenarios in {args.scenarios}")

    want = args.only or list(scenarios)
    unknown = [w for w in want if w not in scenarios]
    if unknown:
        die(f"unknown scenario id(s): {', '.join(unknown)}")

    avail = free_gib(args.out_dir)
    log(f"free space at {args.out_dir}: {avail:.0f} GiB")
    if avail < args.min_free_gib and not args.dry_run:
        die(f"only {avail:.0f} GiB free, want >= {args.min_free_gib:.0f}. The full set is "
            f"~45 GB plus whatever the host run needs.")

    rows: list[dict] = []
    for sid in want:
        sc = scenarios[sid]
        accs = sc.get("accessions") or {}
        if not accs:
            log(f"SKIP {sid}: no `accessions` block (it is a mixture recipe, so fetch the "
                f"scenarios it draws from instead)")
            continue
        declared = sc.get("bracken_readlen")
        truth = sc.get("truth") or {}
        present = ";".join(truth.get("present") or [])
        absent = ";".join(truth.get("absent") or [])
        if not present:
            die(f"{sid} declares no truth.present. A manifest row with empty truth scores "
                f"every call as 'unlisted' and silently measures nothing.")

        for platform, acc_list in accs.items():
            for acc in acc_list or []:
                # ONT read length is not fixed, so the Bracken check does not apply -
                # and Bracken must not be used for ONT at all (build_kraken2_db.sh:526).
                dl = declared if platform == "illumina" else None
                got = fetch_accession(acc, args.out_dir, dl, args.dry_run)
                rows.append({
                    "sample": f"{sid}_{platform}_{acc}",
                    "scenario_id": sid,
                    "kind": sc.get("kind", "mock_community"),
                    "platform": platform,
                    "role": "sample",
                    "replicate": 1,
                    "accession": acc,
                    "run_title": got["row"].get("sample_title", ""),
                    "instrument": got["row"].get("instrument_model", ""),
                    "read_count": got["row"].get("read_count", ""),
                    "bracken_readlen": declared if platform == "illumina" else "",
                    "truth_present": present,
                    "truth_absent": absent,
                    "truth_vaccine": "",
                    "fastq_1": str(got["files"][0]) if got["files"] else "",
                    "fastq_2": str(got["files"][1]) if len(got["files"]) > 1 else "",
                })

    if not rows:
        die("nothing fetched, so no manifest written. Refusing to emit an empty manifest, "
            "which the evaluator would read as a benchmark that measured nothing.")

    if args.dry_run:
        log(f"dry run complete - {len(rows)} manifest row(s) would be written")
        return 0

    cols = list(rows[0])
    man = args.out_dir / "manifest.tsv"
    man.parent.mkdir(parents=True, exist_ok=True)
    with man.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    log(f"wrote {man} ({len(rows)} row(s))")

    sheet = args.out_dir / "samplesheet.csv"
    with sheet.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "platform", "fastq_1", "fastq_2", "role"])
        for r in rows:
            w.writerow([r["sample"], r["platform"], r["fastq_1"], r["fastq_2"], r["role"]])
    log(f"wrote {sheet}")
    log("NOTE: per-run bracken_readlen is in the manifest but the pipeline takes ONE "
        "--bracken_readlen per run, so group samples by read length or run them separately.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
