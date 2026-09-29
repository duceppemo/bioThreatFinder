#!/usr/bin/env python3
"""
Regression tests for bioThreatFinder's decision logic.

WHY THIS FILE EXISTS

Every behaviour asserted here was originally verified once, by hand, in a scratch
directory, and then described in a commit message. That is not tracking validation - it
is remembering it. A claim like "anthrax is held back when plasmid markers are missing"
is only true until someone edits a threshold, and a commit message cannot notice.

So each case below is one behaviour the pipeline is *supposed* to have, with the inputs
that prove it. They are fixture-driven and need no databases, so they run in seconds and
can gate a commit.

  tests/test_pipeline.py            run everything
  tests/test_pipeline.py --list     show the cases
  tests/test_pipeline.py -k marker  run a subset

What this does NOT do: measure accuracy. No ground-truth data is involved. That is
docs/benchmarking.md's job and it has not been run yet. These tests check that the logic
behaves as designed, not that the design is right.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / "bin"
SCRIPTS = REPO / "scripts"
PACK = REPO / "packs" / "human"

COV_COLS = [
    # "platform" is here because score_evidence.py selects the read-identity floor by
    # platform as well as group. It was missing, and cov_row() silently dropped unknown
    # keys, so a fixture passing platform= exercised nothing.
    "sample", "role", "taxid", "taxon", "group", "platform",
    "reference_bp", "reference_contigs",
    "unique_reads", "reads_per_million", "aligned_bases", "breadth", "mean_depth",
    "depth_cv", "window_bp", "windows_expected", "windows_occupied", "window_occupancy",
    "mean_identity", "best_neighbour_taxid", "best_neighbour_reads", "neighbour_ratio",
    "screen_assigned_reads", "mapped_fraction_of_assigned", "total_reads_screened",
    "mapq_min",
]

results: list[tuple[str, str, str]] = []   # (case, PASS/FAIL/SKIP, detail)


def record(case: str, ok: bool, detail: str = "") -> None:
    results.append((case, "PASS" if ok else "FAIL", detail))


def skip(case: str, why: str) -> None:
    results.append((case, "SKIP", why))


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kw)


def cov_row(**kw) -> str:
    d = {c: "" for c in COV_COLS}
    # Reject unknown keys rather than dropping them. The join only walks COV_COLS, so a
    # misspelled or missing column name used to vanish without trace - a fixture would
    # then assert a behaviour it had not actually set up, and pass.
    unknown = sorted(set(kw) - set(COV_COLS))
    if unknown:
        raise KeyError(f"cov_row() got column(s) not in COV_COLS: {unknown}")
    d.update(kw)
    return "\t".join(str(d[c]) for c in COV_COLS)


def write_cov(path: Path, rows: list[str]) -> None:
    path.write_text("\t".join(COV_COLS) + "\n" + "\n".join(rows) + "\n")


def score(tmp: Path, cov: Path, *, markers: Path | None = None,
          candidates: Path | None = None, missing: Path | None = None,
          extra: list[str] | None = None) -> dict[str, dict]:
    """Run score_evidence.py and return findings keyed by sample."""
    out = tmp / f"findings_{cov.stem}.tsv"
    cmd = [sys.executable, BIN / "score_evidence.py", "--coverage", cov,
           "--watchlist", PACK / "watchlist.resolved.tsv",
           "--pack-yaml", PACK / "pack.yaml", "--out", out]
    if markers:
        cmd += ["--marker-evidence", markers, "--marker-rules", PACK / "marker_rules.yaml"]
    if candidates:
        cmd += ["--candidates", candidates]
    if missing:
        cmd += ["--missing-refs", missing]
    cmd += extra or []
    p = run(cmd)
    if not out.exists():
        raise RuntimeError(f"score_evidence produced nothing:\n{p.stderr[-1500:]}")
    # Keyed by sample AND taxon on purpose: a no_reference row carries the same sample
    # name as the real finding, so keying on sample alone silently returns whichever
    # row happened to come last.
    with out.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out_map: dict = {}
    for r in rows:
        out_map[(r["sample"], r["taxon"])] = r
        out_map.setdefault(r["sample"], r) if r["tier"] != "no_reference" else None
    return out_map


# ---------------------------------------------------------------------------
# 1. Self-tests that already live inside the tools
# ---------------------------------------------------------------------------

def test_selftests(tmp: Path) -> None:
    for name, cmd in [
        ("selftest:btf_yaml", [sys.executable, BIN / "btf_yaml.py", "--self-test"]),
        ("selftest:marker_rules", [sys.executable, BIN / "marker_rules.py", "--self-test"]),
        ("selftest:evaluate_benchmark",
         [sys.executable, BIN / "evaluate_benchmark.py", "--self-test",
          "--scenarios", "x", "--manifest", "x", "--findings", "x", "--out-dir", "x"]),
        ("selftest:build_marker_db", [sys.executable, SCRIPTS / "build_marker_db.py", "--self-test"]),
    ]:
        p = run(cmd)
        record(name, p.returncode == 0, p.stderr.strip().splitlines()[-1] if p.returncode else "")


# ---------------------------------------------------------------------------
# 2. Evidence tiering — the six decision paths
# ---------------------------------------------------------------------------

def test_tiering(tmp: Path) -> None:
    cov = tmp / "tiering.tsv"
    write_cov(cov, [
        # Real: spread across the genome, decisively beats its neighbours.
        cov_row(sample="TP", role="sample", taxid="632", taxon="Yersinia pestis",
                group="bacteria", reference_bp=4650000, unique_reads=820,
                reads_per_million=164.0, breadth=0.3400, mean_depth=0.88, depth_cv=1.15,
                window_bp=1000, windows_expected=4650, windows_occupied=2790,
                window_occupancy=0.6000, mean_identity=0.9950, best_neighbour_taxid="633",
                best_neighbour_reads=12, neighbour_ratio="68.33",
                screen_assigned_reads=1000, mapped_fraction_of_assigned="0.8200",
                total_reads_screened=5000000, mapq_min=10),
        # Artifact: piled on a conserved locus, and the neighbour wins.
        # Reads must clear the per-taxon floor (100 for high_consequence bacteria after
        # the UKHSA-derived recalibration) so the case actually reaches the
        # neighbour-attribution gate rather than stopping at the read floor.
        cov_row(sample="ARTIFACT", role="sample", taxid="1392", taxon="Bacillus anthracis",
                group="bacteria", reference_bp=5500000, unique_reads=150,
                reads_per_million=30.0, breadth=0.0021, mean_depth=0.01, depth_cv=7.80,
                window_bp=1000, windows_expected=5500, windows_occupied=9,
                window_occupancy=0.0016, mean_identity=0.9910, best_neighbour_taxid="1396",
                best_neighbour_reads=380, neighbour_ratio="0.3947",
                screen_assigned_reads=400, mapped_fraction_of_assigned="0.3750",
                total_reads_screened=5000000, mapq_min=10),
        # Below the per-taxon read floor.
        cov_row(sample="WEAK", role="sample", taxid="5833", taxon="Plasmodium falciparum",
                group="parasite", reference_bp=23300000, unique_reads=11,
                reads_per_million=2.2, breadth=0.0001, mean_depth=0.0, depth_cv=9.10,
                window_bp=1000, windows_expected=23300, windows_occupied=3,
                window_occupancy=0.0001, mean_identity=0.9300, best_neighbour_taxid="5855",
                best_neighbour_reads=2, neighbour_ratio="5.5",
                screen_assigned_reads=14, mapped_fraction_of_assigned="0.7857",
                total_reads_screened=5000000, mapq_min=10),
        # Present in the negative control at a comparable level.
        cov_row(sample="SAMPLE_CONTAM", role="sample", taxid="2697049",
                taxon="Severe acute respiratory syndrome coronavirus 2", group="virus",
                reference_bp=29900, unique_reads=140, reads_per_million=28.0,
                breadth=0.4200, mean_depth=0.63, depth_cv=1.40, window_bp=598,
                windows_expected=50, windows_occupied=26, window_occupancy=0.5200,
                mean_identity=0.9970, neighbour_ratio="inf",
                screen_assigned_reads=150, mapped_fraction_of_assigned="0.9333",
                total_reads_screened=5000000, mapq_min=10),
        cov_row(sample="NTC", role="control", taxid="2697049",
                taxon="Severe acute respiratory syndrome coronavirus 2", group="virus",
                reference_bp=29900, unique_reads=55, reads_per_million=22.0,
                breadth=0.1900, mean_depth=0.23, depth_cv=1.60, window_bp=598,
                windows_expected=50, windows_occupied=11, window_occupancy=0.2200,
                mean_identity=0.9960, neighbour_ratio="inf",
                screen_assigned_reads=60, mapped_fraction_of_assigned="0.9166",
                total_reads_screened=2500000, mapq_min=10),
    ])
    miss = tmp / "missing.tsv"
    miss.write_text("TP\t10255\tVariola virus\tno_reference\n")
    f = score(tmp, cov, missing=miss)

    expect = {
        "TP": "confirmed_in_silico",
        "ARTIFACT": "attributed_to_neighbour",
        "WEAK": "screened_only",
        "SAMPLE_CONTAM": "suppressed_control",
        "NTC": "control_background",
    }
    for sample, want in expect.items():
        got = f.get(sample, {}).get("tier", "<missing>")
        record(f"tier:{sample}", got == want, f"got {got}, want {want}")

    # ---- the identity floor must be platform-aware ----
    #
    # score_evidence.py had no platform term, so the pack's 0.98 Illumina floor was
    # applied to ONT reads, whose mean identity is ~0.95-0.98. Every bacterial and
    # fungal ONT finding failed on sequencing error alone and capped at `supported`,
    # with a reason string indistinguishable from genuinely weak evidence. Same
    # evidence, both platforms, at an identity ONT really delivers.
    cov_p = tmp / "platform.tsv"
    def pestis(sample: str, platform: str) -> str:
        return cov_row(
            sample=sample, role="sample", taxid="632", taxon="Yersinia pestis",
            group="bacteria", platform=platform, reference_bp=4653728,
            unique_reads=4000, reads_per_million=2000.0, breadth=0.4000,
            mean_depth=5.0, depth_cv=1.00, window_bp=1000, windows_expected=4654,
            windows_occupied=3723, window_occupancy=0.8000,
            mean_identity=0.9520, neighbour_ratio="20.0",
            screen_assigned_reads=4400, mapped_fraction_of_assigned="0.9090",
            total_reads_screened=5000000, mapq_min=10)
    write_cov(cov_p, [pestis("ONT", "ont"), pestis("ILMN", "illumina")])
    fp = score(tmp, cov_p)
    record("tier:ONT bacteria can reach confirmed at ONT-realistic identity",
           fp.get("ONT", {}).get("tier") == "confirmed_in_silico",
           f"tier={fp.get('ONT', {}).get('tier')} "
           f"reasons={fp.get('ONT', {}).get('reasons', '')[:120]}")
    record("tier:Illumina still holds the strict floor at the same identity",
           fp.get("ILMN", {}).get("tier") == "supported",
           f"tier={fp.get('ILMN', {}).get('tier')}")
    record("tier:the reason names which platform rule was applied",
           "ont/default" in fp.get("ONT", {}).get("reasons", ""),
           fp.get("ONT", {}).get("reasons", "")[:160])

    # A database gap must NOT be reported as a negative.
    gap = [r for r in f.values() if r["tier"] == "no_reference"]
    record("tier:no_reference surfaces as a database gap", len(gap) == 1,
           f"{len(gap)} no_reference rows")


# ---------------------------------------------------------------------------
# 3. Marker gating — identical alignment evidence, four outcomes
# ---------------------------------------------------------------------------

def test_marker_gating(tmp: Path) -> None:
    def anthrax(sample: str) -> str:
        # Alignment evidence that WOULD confirm on its own, held constant.
        return cov_row(sample=sample, role="sample", taxid="1392",
                       taxon="Bacillus anthracis", group="bacteria",
                       reference_bp=5500000, unique_reads=900, reads_per_million=180.0,
                       breadth=0.3800, mean_depth=0.95, depth_cv=1.10, window_bp=1000,
                       windows_expected=5500, windows_occupied=3300,
                       window_occupancy=0.6000, mean_identity=0.9900,
                       best_neighbour_taxid="1396", best_neighbour_reads=14,
                       neighbour_ratio="64.29", screen_assigned_reads=1000,
                       mapped_fraction_of_assigned="0.9000",
                       total_reads_screened=5000000, mapq_min=10)

    cov = tmp / "gating.tsv"
    write_cov(cov, [
        anthrax("BOTH"), anthrax("PXO1_ONLY"), anthrax("UNTESTED"), anthrax("VETOED"),
        # MTBC with RD1 deleted = BCG vaccine strain, not disease.
        cov_row(sample="BCG", role="sample", taxid="1773",
                taxon="Mycobacterium tuberculosis", group="bacteria",
                reference_bp=4410000, unique_reads=600, reads_per_million=120.0,
                breadth=0.3500, mean_depth=0.82, depth_cv=1.20, window_bp=1000,
                windows_expected=4410, windows_occupied=2600, window_occupancy=0.5900,
                mean_identity=0.9920, best_neighbour_taxid="1765", best_neighbour_reads=9,
                neighbour_ratio="66.67", screen_assigned_reads=700,
                mapped_fraction_of_assigned="0.8571", total_reads_screened=5000000,
                mapq_min=10),
    ])
    mk = tmp / "markers.tsv"
    mk.write_text("sample\tmarker_id\tstate\n" + "\n".join([
        "BOTH\tbanth_pxo1_pagA\tPRESENT", "BOTH\tbanth_pxo2_capB\tPRESENT",
        "BOTH\tbcereus_chr_ces\tABSENT",
        "PXO1_ONLY\tbanth_pxo1_pagA\tPRESENT", "PXO1_ONLY\tbanth_pxo2_capB\tABSENT",
        "PXO1_ONLY\tbanth_pxo2_capA\tABSENT", "PXO1_ONLY\tbcereus_chr_ces\tABSENT",
        # UNTESTED: pXO2 markers deliberately ABSENT FROM THE FILE, not ABSENT in value.
        "UNTESTED\tbanth_pxo1_pagA\tPRESENT",
        "VETOED\tbanth_pxo1_pagA\tPRESENT", "VETOED\tbanth_pxo2_capB\tPRESENT",
        "VETOED\tbcereus_chr_ces\tPRESENT",
        "BCG\tmtb_rpoB\tPRESENT", "BCG\tmtb_RD9\tABSENT", "BCG\tmtb_RD1\tABSENT",
    ]) + "\n")
    f = score(tmp, cov, markers=mk)

    for sample, want in {
        "BOTH": "confirmed_in_silico",
        "PXO1_ONLY": "supported",       # required markers confidently absent
        "UNTESTED": "supported",        # untested != absent
        "VETOED": "marker_veto",
        "BCG": "vaccine_strain",
    }.items():
        got = f.get(sample, {}).get("tier", "<missing>")
        record(f"marker:{sample}", got == want, f"got {got}, want {want}")

    # The untested case must SAY it is a database gap, not just score low.
    why = f.get("UNTESTED", {}).get("reasons", "")
    record("marker:untested names the gap explicitly",
           "INCOMPLETE" in why and "database gap" in why.lower(), why[:120])
    # A confirmation whose veto could not be evaluated must disclose that.
    cov2 = tmp / "gating2.tsv"
    write_cov(cov2, [anthrax("NOVETO")])
    mk2 = tmp / "markers2.tsv"
    mk2.write_text("sample\tmarker_id\tstate\nNOVETO\tbanth_pxo1_pagA\tPRESENT\n"
                   "NOVETO\tbanth_pxo2_capB\tPRESENT\n")
    f2 = score(tmp, cov2, markers=mk2)
    record("marker:unevaluated veto is disclosed",
           "CAVEAT" in f2.get("NOVETO", {}).get("reasons", ""),
           f2.get("NOVETO", {}).get("tier", ""))
    # AND it must cap the tier, not merely mention it in prose. B. anthracis's veto
    # marker (bcereus_chr_ces) is untested here, so the near-neighbour argument against
    # the call was never made - confirming anyway treats an ungathered check as passed.
    record("marker:an unevaluable veto HOLDS the finding at supported",
           f2.get("NOVETO", {}).get("tier") == "supported",
           f"tier={f2.get('NOVETO', {}).get('tier')}")

    # The vaccine case, which is the one with a public-health consequence. 18 markers
    # referenced by rules can only ever be UNKNOWN, so 8 of 9 vaccine_if rules can never
    # fire; a Sabin vaccine strain satisfied Enterovirus C's confirm_requires while its
    # polio_sabin check was untestable, and was reported as confirmed poliovirus.
    cov3 = tmp / "gating3.tsv"
    write_cov(cov3, [cov_row(
        sample="SABIN", role="sample", taxid="138950", taxon="Enterovirus C",
        group="virus", reference_bp=7440, unique_reads=900, reads_per_million=180.0,
        breadth=0.9500, mean_depth=12.0, depth_cv=0.60, window_bp=1000,
        windows_expected=8, windows_occupied=8, window_occupancy=1.0000,
        mean_identity=0.9900, best_neighbour_taxid="12080",
        best_neighbour_reads=5, neighbour_ratio="180.0",
        screen_assigned_reads=1000, mapped_fraction_of_assigned="0.9000",
        total_reads_screened=5000000, mapq_min=10, platform="illumina")]) 
    mk3 = tmp / "markers3.tsv"
    # polio_VP1 present (built), polio_sabin absent from the FILE = untested.
    mk3.write_text("sample\tmarker_id\tstate\nSABIN\tpolio_VP1\tPRESENT\n")
    f3 = score(tmp, cov3, markers=mk3)
    sabin_tier = f3.get("SABIN", {}).get("tier", "<missing>")
    sabin_why = f3.get("SABIN", {}).get("reasons", "")
    record("marker:an untestable vaccine check HOLDS the finding at supported",
           sabin_tier != "confirmed_in_silico",
           f"tier={sabin_tier} - confirming here reports a vaccine strain as "
           f"wild-type disease")
    record("marker:and says a vaccine strain could not be excluded",
           "vaccine strain CANNOT be excluded" in sabin_why, sabin_why[:160])


# ---------------------------------------------------------------------------
# 4. Wright et al. confirmation composite
# ---------------------------------------------------------------------------

def test_wright_gates(tmp: Path) -> None:
    cov = tmp / "wright.tsv"
    write_cov(cov, [
        cov_row(sample="GENUINE", role="sample", taxid="632", taxon="Yersinia pestis",
                group="bacteria", reference_bp=4650000, unique_reads=820,
                reads_per_million=164.0, breadth=0.3400, mean_depth=0.88, depth_cv=1.15,
                window_bp=1000, windows_expected=4650, windows_occupied=2790,
                window_occupancy=0.6000, mean_identity=0.9950, best_neighbour_taxid="633",
                best_neighbour_reads=12, neighbour_ratio="68.33",
                screen_assigned_reads=1000, mapped_fraction_of_assigned="0.8200",
                total_reads_screened=5000000, mapq_min=10),
        # Wright's Elizabethkingia signature: even coverage, good breadth, ~96% identity.
        # Genus right, SPECIES WRONG. The old 0.90 floor confirmed this.
        cov_row(sample="NEAR_NEIGHBOUR", role="sample", taxid="446",
                taxon="Legionella pneumophila", group="bacteria", reference_bp=3400000,
                unique_reads=600, reads_per_million=120.0, breadth=0.3000,
                mean_depth=0.70, depth_cv=1.20, window_bp=1000, windows_expected=3400,
                windows_occupied=1870, window_occupancy=0.5500, mean_identity=0.9620,
                best_neighbour_taxid="452", best_neighbour_reads=30,
                neighbour_ratio="20.00", screen_assigned_reads=660,
                mapped_fraction_of_assigned="0.9090", total_reads_screened=5000000,
                mapq_min=10),
        # k-mer evidence did not survive alignment: only 18% of assigned reads map.
        cov_row(sample="LOW_MAPPED", role="sample", taxid="1639",
                taxon="Listeria monocytogenes", group="bacteria", reference_bp=2940000,
                unique_reads=180, reads_per_million=36.0, breadth=0.2500,
                mean_depth=0.30, depth_cv=1.30, window_bp=1000, windows_expected=2940,
                windows_occupied=1470, window_occupancy=0.5000, mean_identity=0.9900,
                best_neighbour_taxid="1642", best_neighbour_reads=9,
                neighbour_ratio="20.00", screen_assigned_reads=1000,
                mapped_fraction_of_assigned="0.1800", total_reads_screened=5000000,
                mapq_min=10),
        # Genuinely divergent virus at 89.2%. A bacterial floor here is a false negative.
        cov_row(sample="DIVERGENT_VIRUS", role="sample", taxid="3052518",
                taxon="Crimean-Congo hemorrhagic fever orthonairovirus", group="virus",
                reference_bp=12000, reference_contigs=3, unique_reads=150,
                reads_per_million=30.0, breadth=0.3500, mean_depth=0.90, depth_cv=1.40,
                window_bp=240, windows_expected=50, windows_occupied=30,
                window_occupancy=0.6000, mean_identity=0.8920, neighbour_ratio="inf",
                screen_assigned_reads=190, mapped_fraction_of_assigned="0.7890",
                total_reads_screened=5000000, mapq_min=10),
    ])
    mk = tmp / "wright_markers.tsv"
    mk.write_text("sample\tmarker_id\tstate\n" + "\n".join([
        "GENUINE\typest_ppcp1_pla\tPRESENT", "GENUINE\typest_pmt1_caf1\tPRESENT",
        "GENUINE\typseudo_chr_inv\tABSENT",
        "NEAR_NEIGHBOUR\tlpneu_mip\tPRESENT",
        "LOW_MAPPED\tlmono_hly\tPRESENT", "LOW_MAPPED\tlmono_prfA\tPRESENT",
        "DIVERGENT_VIRUS\tcchf_L\tPRESENT",
    ]) + "\n")
    f = score(tmp, cov, markers=mk)
    for sample, want in {
        "GENUINE": "confirmed_in_silico",
        "NEAR_NEIGHBOUR": "supported",
        "LOW_MAPPED": "supported",
        "DIVERGENT_VIRUS": "confirmed_in_silico",
    }.items():
        got = f.get(sample, {}).get("tier", "<missing>")
        record(f"wright:{sample}", got == want, f"got {got}, want {want}")

    record("wright:96% identity fails on the identity gate",
           "identity" in f.get("NEAR_NEIGHBOUR", {}).get("reasons", ""),
           f.get("NEAR_NEIGHBOUR", {}).get("reasons", "")[:110])
    record("wright:low mapped fraction fails on that gate",
           "mapped_fraction" in f.get("LOW_MAPPED", {}).get("reasons", ""),
           f.get("LOW_MAPPED", {}).get("reasons", "")[:110])


# ---------------------------------------------------------------------------
# 5. Coverage metrics against a hand-built BAM, incl. masking
# ---------------------------------------------------------------------------

def test_coverage_metrics(tmp: Path) -> None:
    if not shutil.which("samtools"):
        skip("coverage:metrics", "samtools not on PATH")
        skip("coverage:masking", "samtools not on PATH")
        return

    tgt = "taxid|632|NZ_TEST1.1 Yersinia pestis chr"
    nb = "taxid|633|NC_TEST2.1 Yersinia pseudotuberculosis"
    L = 1000

    def rd(name, ref, pos, mapq, flag=0, nm=1, ln=100):
        return (f"{name}\t{flag}\t{ref}\t{pos}\t{mapq}\t{ln}M\t*\t0\t0"
                f"\t{'A'*ln}\t{'I'*ln}\tNM:i:{nm}")

    sam = tmp / "t.sam"
    lines = ["@HD\tVN:1.6\tSO:coordinate", f"@SQ\tSN:{tgt}\tLN:{L}", f"@SQ\tSN:{nb}\tLN:{L}"]
    lines += [rd(f"f{i}", tgt, p, 60) for i, p in enumerate([1, 201, 401, 601])]
    lines += [rd("r0", tgt, 801, 60, flag=16),      # reverse strand
              rd("dup0", tgt, 1, 60),                # duplicate start position
              rd("lowq", tgt, 301, 5),               # below the MAPQ floor
              rd("nb0", nb, 1, 60)]                  # the competing neighbour
    sam.write_text("\n".join(lines) + "\n")
    bam = tmp / "t.bam"
    if run(["samtools", "sort", "-o", bam, sam]).returncode != 0:
        skip("coverage:metrics", "samtools sort failed")
        return
    run(["samtools", "index", bam])

    def cov_of(extra=None):
        out = tmp / f"cov{'_masked' if extra else ''}.tsv"
        p = run([sys.executable, BIN / "coverage_stats.py", "--bam", bam, "--sample", "T",
                 "--taxid", "632", "--taxon", "Yersinia pestis", "--group", "bacteria",
                 "--mapq-min", "10", "--total-reads", "1000000", "--assigned-reads", "10",
                 "--out", out] + (extra or []))
        if not out.exists():
            raise RuntimeError(p.stderr[-800:])
        with out.open() as fh:
            return next(csv.DictReader(fh, delimiter="\t"))

    r = cov_of()
    for key, want in {
        "unique_reads": "6",                    # 4 fwd + 1 rev + 1 dup; lowq excluded
        "mapped_reads_any_mapq": "7",           # + the low-MAPQ read
        "distinct_start_positions": "5",
        "forward_proportion": "0.8333",         # 5 of 6
        "duplication_rate_proxy": "0.1667",
        "mapq_pass_fraction": "0.8571",
        "best_neighbour_taxid": "633",
        "mapped_fraction_of_assigned": "0.6000",
        "mean_alignment_length": "100.0",
        "breadth": "0.500000",
    }.items():
        record(f"coverage:{key}", r.get(key) == want, f"got {r.get(key)}, want {want}")

    # Masking: 140 bp masked, fully containing the read at position 601.
    mdir = tmp / "mask"
    mdir.mkdir()
    (mdir / "t.mask.bed").write_text("NZ_TEST1.1\t580\t720\tcrossgroup\n")
    m = cov_of(["--mask-bed", mdir])
    for key, want in {
        "unique_reads": "5",                    # the masked read no longer counts
        "reads_in_masked_regions": "1",
        "masked_bp": "140",
        "reference_bp": "860",                  # denominator shrinks
        "reference_bp_raw": "1000",
        "breadth": "0.465116",
    }.items():
        record(f"coverage:masked:{key}", m.get(key) == want, f"got {m.get(key)}, want {want}")

    # ---- occupancy must be a FRACTION on a segmented reference ----
    #
    # window_hits is keyed per contig over raw coordinates, while the denominator was
    # ceil(mask-subtracted total / window) across ALL contigs. Each contig contributes
    # ceil(len_i/window) windows and their sum exceeds the pooled ceiling, so a fully
    # covered 8-segment influenza-shaped reference scored occupancy 1.12 - loosening
    # the occupancy gate in the false-positive direction on exactly the segmented
    # viruses and fragmented drafts where it should bite hardest.
    seg_len, n_seg = 1400, 8
    seg_refs = [f"taxid|11320|SEG{i}.1 Influenza A segment {i}" for i in range(n_seg)]
    slines = ["@HD\tVN:1.6\tSO:coordinate"]
    slines += [f"@SQ\tSN:{r}\tLN:{seg_len}" for r in seg_refs]
    # One full-length read per segment: every window of every segment is covered.
    slines += [rd(f"seg{i}", r, 1, 60, ln=seg_len) for i, r in enumerate(seg_refs)]
    ssam = tmp / "seg.sam"
    ssam.write_text("\n".join(slines) + "\n")
    sbam = tmp / "seg.bam"
    if run(["samtools", "sort", "-o", sbam, ssam]).returncode == 0:
        run(["samtools", "index", sbam])
        sout = tmp / "seg.tsv"
        run([sys.executable, BIN / "coverage_stats.py", "--bam", sbam, "--sample", "T",
             "--taxid", "11320", "--taxon", "Influenza A virus", "--group", "virus",
             "--platform", "illumina", "--mapq-min", "10", "--total-reads", "1000000",
             "--assigned-reads", "8", "--out", sout])
        if sout.exists():
            with sout.open() as fh:
                sr = next(csv.DictReader(fh, delimiter="\t"))
            occ = float(sr.get("window_occupancy") or 0)
            record("coverage:occupancy never exceeds 1.0 on a segmented reference",
                   occ <= 1.0, f"occupancy={occ} across {n_seg} segments")
            record("coverage:a fully covered segmented reference is fully occupied",
                   abs(occ - 1.0) < 1e-6, f"occupancy={occ}, want 1.0")
            record("coverage:platform is recorded for the identity floor",
                   sr.get("platform") == "illumina", f"got {sr.get('platform')!r}")


# ---------------------------------------------------------------------------
# 6. Multi-classifier screen parsing and the unclassified fraction
# ---------------------------------------------------------------------------

def test_screen(tmp: Path) -> None:
    desc = tmp / "desc.tsv"
    p = run([sys.executable, SCRIPTS / "resolve_watchlist.py",
             "--watchlist", PACK / "watchlist.tsv",
             "--taxdump", os.environ.get("BTF_TAXDUMP", "/media/8tb/bioThreatFinder/taxdump/2026-09-29"),
             "--out", tmp / "wl.tsv", "--descendants-out", desc])
    if p.returncode != 0 or not desc.exists():
        skip("screen:multi-classifier", "no usable taxdump (set BTF_TAXDUMP)")
        skip("screen:unclassified fraction", "no usable taxdump")
        return
    wl = tmp / "wl.tsv"

    # KrakenUniq: identical read counts, opposite distinct-k-mer evidence.
    ku = tmp / "ku.tsv"
    ku.write_text("\n".join([
        "%\treads\ttaxReads\tkmers\tdup\tcov\ttaxID\trank\ttaxName",
        "0.004\t800\t760\t45000\t1.18\t0.312\t632\tspecies\tYersinia pestis",
        "0.003\t600\t580\t78\t41.2\t0.0009\t1392\tspecies\tBacillus anthracis",
    ]) + "\n")
    # Kraken2 with the unclassified + root rows, 55% unclassified.
    k2 = tmp / "k2.tsv"
    k2.write_text("\n".join([
        "55.0000\t550000\t550000\tU\t0\tunclassified",
        "45.0000\t450000\t0\tR\t1\troot",
        "0.0800\t800\t760\tS\t632\t  Yersinia pestis",
    ]) + "\n")
    cand, summ = tmp / "cand.tsv", tmp / "summary.tsv"
    p = run([sys.executable, BIN / "parse_screen.py", "--krakenuniq-report", ku,
             "--kraken2-report", k2, "--watchlist", wl, "--descendants", desc,
             "--sample", "S1", "--out", cand, "--out-summary", summ])
    if not cand.exists():
        record("screen:multi-classifier", False, p.stderr[-200:])
        return
    with cand.open() as fh:
        rows = {r["name"]: r for r in csv.DictReader(fh, delimiter="\t")}

    record("screen:clean taxon carries high distinct k-mers",
           rows.get("Yersinia pestis", {}).get("krakenuniq_distinct_kmers") == "45000")
    anth = rows.get("Bacillus anthracis", {})
    record("screen:piled-up taxon is flagged, not dropped",
           "low_distinct_kmers" in anth.get("screen_flags", "")
           and "high_kmer_duplication" in anth.get("screen_flags", ""),
           anth.get("screen_flags", "<no row>"))
    record("screen:HLL cardinality read from the clade, never summed",
           anth.get("kmer_basis") == "clade", anth.get("kmer_basis", ""))

    if summ.exists():
        with summ.open() as fh:
            s = {r["classifier"]: r for r in csv.DictReader(fh, delimiter="\t")}
        record("screen:unclassified fraction is reported and flagged",
               s.get("kraken2", {}).get("flag") == "HIGH_UNCLASSIFIED"
               and s.get("kraken2", {}).get("unclassified_fraction", "").startswith("0.55"),
               str(s.get("kraken2", {})))
    else:
        record("screen:unclassified fraction is reported and flagged", False, "no summary")


# ---------------------------------------------------------------------------
# 7. Pack integrity — the things whose breakage is silent
# ---------------------------------------------------------------------------

def test_pack_integrity(tmp: Path) -> None:
    p = run([sys.executable, BIN / "marker_rules.py", "--rules", PACK / "marker_rules.yaml",
             "--markers", PACK / "markers.tsv", "--validate"])
    record("pack:rules compile and reference only real markers", p.returncode == 0,
           p.stderr.strip().splitlines()[-1] if p.returncode else "")

    sys.path.insert(0, str(BIN))
    try:
        import btf_yaml
        doc = btf_yaml.load(REPO / "benchmarks" / "scenarios" / "human_core.yaml")
        sc = doc.get("scenarios") or []
        record("pack:scenario file parses", len(sc) > 0, f"{len(sc)} scenarios")
        need = {"dna_virus", "rna_virus", "gram_positive_bacteria",
                "gram_negative_bacteria", "mold", "yeast", "parasite"}
        have = {s.get("nysdoh_category") for s in sc}
        record("pack:all 7 NYSDOH pathogen categories covered", need <= have,
               f"missing {sorted(need - have)}")
        record("pack:a plumbing positive control exists",
               any(s.get("id") == "plumbing_positive_control" for s in sc))
        # Every taxon named in a scenario must have a genome source in the pack.
        with (PACK / "watchlist.tsv").open() as fh:
            wl_rows = [r for r in csv.DictReader(
                [l for l in fh if not l.startswith("#")], delimiter="\t") if r.get("name")]
        known = {r["name"] for r in wl_rows}
        for r in wl_rows:
            known |= {x.strip() for x in (r.get("near_neighbours") or "").split(";")
                      if x.strip() and x.strip() != "-"}
        named = set((doc.get("defaults") or {}).get("background") or {})
        for s in sc:
            for sp in s.get("spikes") or []:
                named.add(sp["taxon"])
            for k in ("present", "absent", "vaccine"):
                named |= set((s.get("truth") or {}).get(k) or [])
        record("pack:every scenario taxon has a genome source", not (named - known),
               f"orphans: {sorted(named - known)}")
        pk = btf_yaml.load(PACK / "pack.yaml")
        th = pk.get("thresholds") or {}
        record("pack:identity floor is the near-neighbour discriminator (>=0.98)",
               float(th.get("min_read_identity", 0)) >= 0.98,
               str(th.get("min_read_identity")))
        record("pack:viruses get a lower identity floor than bacteria",
               float((th.get("min_read_identity_by_group") or {}).get("virus", 1.0))
               < float(th.get("min_read_identity", 0)))

        # A THRESHOLD MUST BE REACHABLE, NOT MERELY STRICT.
        #
        # The two checks above assert the floor is high, which is what let the ONT
        # defect hide: 0.98 applied to ONT reads, whose mean identity is ~0.95-0.98,
        # made every bacterial and fungal ONT confirmation impossible while the suite
        # stayed green. These assert the floor can actually be met on each platform.
        import score_evidence as se
        # badread's default simulated identity, i.e. what the ONT benchmark delivers.
        ONT_SIMULATED_MEAN = 0.95
        for grp in ("bacteria", "fungus", "virus", "parasite"):
            floor, src = se.resolve_identity_floor(th, grp, "ont")
            record(f"pack:ONT {grp} identity floor is reachable (< {ONT_SIMULATED_MEAN})",
                   floor < ONT_SIMULATED_MEAN,
                   f"floor={floor} via {src}; ONT reads cannot exceed ~"
                   f"{ONT_SIMULATED_MEAN} as simulated")
        il, il_src = se.resolve_identity_floor(th, "bacteria", "illumina")
        record("pack:Illumina keeps the strict near-neighbour floor (>=0.98)",
               il >= 0.98, f"floor={il} via {il_src}")
        ont_b, _ = se.resolve_identity_floor(th, "bacteria", "ont")
        record("pack:ONT floor is looser than Illumina for the same group",
               ont_b < il, f"ont={ont_b} illumina={il}")
        _, nosrc = se.resolve_identity_floor(th, "bacteria", "")
        record("pack:a missing platform is named in the reason string",
               "no platform" in nosrc, nosrc)
    except Exception as exc:
        record("pack:scenario file parses", False, f"{type(exc).__name__}: {exc}")


# ---------------------------------------------------------------------------
# 8. Benchmark evaluator: the assertions that gate CI must be able to fail
# ---------------------------------------------------------------------------

def test_evaluator(tmp: Path) -> None:
    sys.path.insert(0, str(BIN))
    try:
        import evaluate_benchmark as eb
    except Exception as exc:
        record("eval:importable", False, f"{type(exc).__name__}: {exc}")
        return

    def ob(tier, truth="present", sid="s1", taxon="T", sample="x"):
        return {"sample": sample, "scenario_id": sid, "taxon": taxon,
                "platform": "illumina", "level_rpm": "100", "truth": truth,
                "tier": tier, "role": "sample", "outcome": eb.outcome_of(tier)}

    def result_of(exp, obs):
        rs = eb.check_expectations({"s1": {"expect": [dict(taxon="T", **exp)]}}, obs)
        return rs[0]["result"] if rs else "NONE"

    # tier_rank() returns -1 for anything unrecognised and nothing ranks below -1, so
    # a misspelled floor could never fail - it removed the gate it was written to be.
    record("eval:a misspelled min_tier FAILS rather than passing silently",
           result_of({"min_tier": "vaccine_strainX"}, [ob("no_reference")]) == "FAIL")
    record("eval:a misspelled max_tier FAILS rather than passing silently",
           result_of({"max_tier": "supportedX"}, [ob("confirmed_in_silico")]) == "FAIL")

    # A missing reference is a database gap, not a threshold regression.
    record("eval:a NOT_ASSESSED observation SKIPs a min_tier check",
           result_of({"min_tier": "supported"}, [ob("no_reference")]) == "SKIP")
    record("eval:a real shortfall still FAILs a min_tier check",
           result_of({"min_tier": "confirmed_in_silico"}, [ob("supported")]) == "FAIL")
    record("eval:a met min_tier still PASSes",
           result_of({"min_tier": "supported"}, [ob("confirmed_in_silico")]) == "PASS")

    # vaccine_table was the one table that counted NOT_ASSESSED in its denominator,
    # so a database gap reported vaccine discrimination as totally broken.
    vt = eb.vaccine_table([ob("no_reference", truth="vaccine")])[0]
    record("eval:vaccine accuracy is blank, not 0.000, when nothing is assessable",
           vt["n_assessable"] == 0 and vt["vaccine_accuracy"] == "",
           f"{vt['n_assessable']=} {vt['vaccine_accuracy']=}")

    # On a taxon the truth says is absent, "vaccine strain of X is present" is a
    # false positive; counting only POSITIVE let it improve the specificity figure.
    sp = eb.specificity_table([ob("vaccine_strain", truth="absent")])[0]
    record("eval:a vaccine call on an absent taxon counts as a false positive",
           sp["n_false_positive"] == 1, f"got {sp['n_false_positive']}")
    vs = eb.variant_stats([ob("vaccine_strain", truth="absent")])
    record("eval:variant_stats counts it too", vs["false_positive"] == 1,
           f"got {vs['false_positive']}")

    # Two scenarios measuring one taxon must not collapse into one row.
    rows = eb.bakeoff_per_taxon(
        {"v1": [ob("supported", sid="a"), ob("supported", sid="b")]})
    sids = {r["scenario_id"] for r in rows}
    record("eval:bake-off rows keep scenario_id so same-taxon scenarios do not collide",
           sids == {"a", "b"}, f"scenario ids in output: {sorted(sids)}")


CASES = {
    "selftests": test_selftests,
    "tiering": test_tiering,
    "marker": test_marker_gating,
    "wright": test_wright_gates,
    "coverage": test_coverage_metrics,
    "screen": test_screen,
    "pack": test_pack_integrity,
    "eval": test_evaluator,
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-k", "--only", help="substring match on the group name")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    if args.list:
        for name in CASES:
            print(f"  {name}")
        return 0

    with tempfile.TemporaryDirectory(prefix="btf_tests.") as td:
        tmp = Path(td)
        for name, fn in CASES.items():
            if args.only and args.only not in name:
                continue
            try:
                fn(tmp)
            except Exception as exc:
                record(f"{name}:<crashed>", False, f"{type(exc).__name__}: {exc}")

    width = max((len(c) for c, _, _ in results), default=10)
    for case, status, detail in results:
        mark = {"PASS": "ok  ", "FAIL": "FAIL", "SKIP": "skip"}[status]
        line = f"  {mark} {case:<{width}}"
        if status != "PASS" and detail:
            line += f"  {detail}"
        print(line)

    n_pass = sum(1 for _, s, _ in results if s == "PASS")
    n_fail = sum(1 for _, s, _ in results if s == "FAIL")
    n_skip = sum(1 for _, s, _ in results if s == "SKIP")
    print(f"\n{n_pass} passed, {n_fail} failed, {n_skip} skipped")
    if n_fail:
        print("\nThese assert DESIGNED BEHAVIOUR, not accuracy. A failure means the logic "
              "changed;\nit does not mean the pipeline got less accurate - nothing here "
              "measures accuracy.")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
