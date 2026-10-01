#!/usr/bin/env python3
"""Fetch assemblies that kraken2-build failed to download, ready for --add-to-library.

WHY THIS EXISTS

`kraken2-build --download-library bacteria` printed 5,413 `gzip: No such file or
directory` errors while concatenating, then reported "library complete" and moved on.
Those assemblies failed to transfer and are silently absent from library.fna - 8% of the
67,307-genome target set. Spot-checking three of them returned HTTP 200, so they exist
and the failures were transient FTP errors over a 22-hour download, not withdrawals.

The detection cost is modest: 24 watchlist taxa lost 5-15% of their assemblies and none
lost all of them, and near-identical assemblies of one species contribute few new
minimizers. The provenance cost is the real one - an index that silently omits 8% of its
declared content cannot be described accurately in a validation record.

WHAT THIS DOES NOT DO

It does not touch the kraken2 database directory. Output goes to a separate staging
directory, so this can run alongside a live `kraken2-build` without the two colliding -
which is the hazard that put an flock on the builders in the first place.

Feed the result to `kraken2-build --add-to-library` once the library downloads have
finished and BEFORE the --build step, then re-run the builder: its stamps skip the
completed libraries and go straight to the build.

Headers are rewritten to `kraken:taxid|<taxid>|<original>`, which kraken2 accepts
directly - necessary because this build used --skip-maps, so there is no
accession2taxid available to resolve plain accessions.
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import sys
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def log(msg: str) -> None:
    sys.stderr.write(f"[fetch] {msg}\n")
    sys.stderr.flush()


def read_summary(path: Path, wanted: set[str]) -> dict[str, tuple[str, str]]:
    """accession -> (taxid, url). NCBI assembly_summary: f6 taxid, f20 ftp_path."""
    out: dict[str, tuple[str, str]] = {}
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 20 or f[0] not in wanted:
                continue
            ftp = f[19].strip()
            if not ftp:
                continue
            base = ftp.rstrip("/").split("/")[-1]
            url = f"{ftp.replace('ftp://', 'https://')}/{base}_genomic.fna.gz"
            out[f[0]] = (f[5].strip(), url)
    return out


def verify_gz(path: Path) -> bool:
    """A truncated transfer that still looks like a file is this project's classic bug."""
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return path.stat().st_size > 0
    except (OSError, EOFError, gzip.BadGzipFile):
        return False


def fetch(acc: str, url: str, dest: Path, retries: int = 4) -> tuple[str, bool, str]:
    if dest.is_file() and verify_gz(dest):
        return acc, True, "cached"
    tmp = dest.with_suffix(dest.suffix + ".part")
    last = ""
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "bioThreatFinder"})
            with urllib.request.urlopen(req, timeout=120) as r, tmp.open("wb") as out:
                shutil.copyfileobj(r, out, length=1 << 20)
            if verify_gz(tmp):
                tmp.replace(dest)
                return acc, True, f"ok (attempt {attempt})"
            last = "downloaded but failed gzip integrity"
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            last = f"{type(exc).__name__}: {exc}"
        tmp.unlink(missing_ok=True)
    return acc, False, last


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accessions", required=True, type=Path,
                    help="one assembly accession per line")
    ap.add_argument("--assembly-summary", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path,
                    help="staging directory; must NOT be inside the kraken2 DB")
    ap.add_argument("--combined", type=Path, default=None,
                    help="write one kraken2-ready FASTA here [<out-dir>/missing.fna]")
    ap.add_argument("--jobs", type=int, default=6,
                    help="concurrent downloads; kept modest to be a good NCBI citizen [6]")
    args = ap.parse_args(argv)

    wanted = {l.strip() for l in args.accessions.read_text().splitlines() if l.strip()}
    if not wanted:
        log("ERROR: no accessions given")
        return 2
    if "pluspf" in str(args.out_dir):
        log(f"ERROR: --out-dir {args.out_dir} is inside the kraken2 database. Two "
            f"kraken2-build processes sharing a directory can interleave partial files "
            f"into library.fna, and a corrupt reference in an index does not error - it "
            f"misclassifies. Stage elsewhere.")
        return 2

    meta = read_summary(args.assembly_summary, wanted)
    log(f"{len(wanted)} accession(s) requested, {len(meta)} resolved to a URL")
    unresolved = sorted(wanted - set(meta))
    if unresolved:
        log(f"WARNING: {len(unresolved)} have no ftp_path in the summary "
            f"(withdrawn?): {', '.join(unresolved[:5])}")

    gz_dir = args.out_dir / "gz"
    gz_dir.mkdir(parents=True, exist_ok=True)
    ok, failed = [], []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = {pool.submit(fetch, a, u, gz_dir / f"{a}.fna.gz"): a
                for a, (_t, u) in meta.items()}
        for done, fut in enumerate(as_completed(futs), start=1):
            acc, good, why = fut.result()
            (ok if good else failed).append(acc)
            if done % 250 == 0 or not good:
                log(f"{done}/{len(meta)}  ok={len(ok)} failed={len(failed)}"
                    + ("" if good else f"  <- {acc}: {why}"))

    log(f"downloaded {len(ok)}, failed {len(failed)}")
    if failed:
        (args.out_dir / "failed.txt").write_text("\n".join(sorted(failed)) + "\n")
        log(f"failures listed in {args.out_dir / 'failed.txt'}")

    combined = args.combined or (args.out_dir / "missing.fna")
    n_seq = 0
    tmp = combined.with_suffix(".fna.part")
    with tmp.open("w", encoding="utf-8", newline="\n") as out:
        for acc in sorted(ok):
            taxid = meta[acc][0]
            with gzip.open(gz_dir / f"{acc}.fna.gz", "rt",
                           encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if line.startswith(">"):
                        n_seq += 1
                        out.write(f">kraken:taxid|{taxid}|{line[1:]}")
                    else:
                        out.write(line)
    tmp.replace(combined)
    log(f"wrote {combined} ({n_seq} sequences from {len(ok)} assemblies, "
        f"{combined.stat().st_size / 1e9:.1f} GB)")
    log("next: kraken2-build --add-to-library "
        f"{combined} --db <DB>   (AFTER the library downloads finish, BEFORE --build)")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
