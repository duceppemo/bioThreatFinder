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
import re
import shutil
import subprocess
import sys
import pathlib
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
    "mean_identity",
    # Absent from this list entirely, so no scoring fixture could set it, so
    # `as_float(row.get("mean_alignment_length"), -1.0)` returned -1 and
    # score_evidence.py SKIPPED the gate in every test ever written. Not "tested at the
    # boundary" - never executed. It was simultaneously unreachable on real 2x100 data,
    # which is how the project's entire top tier stayed broken while 110 tests passed.
    "mean_alignment_length",
    # Per-reference breadth, and which reference it was. Only differs from `breadth`
    # on a bundle holding more than one organism, which is what a group-rank entry is.
    "best_reference_breadth", "best_reference_occupancy", "best_reference",
    "best_neighbour_taxid", "best_neighbour_reads",
    # The D12 credibility gate's ONLY input. It was absent from this list while
    # being present in every real coverage TSV, so cov_row() rejected it and no
    # fixture could reach either side of the gate - the branch that withholds an
    # attribution was unreachable from the test suite that is supposed to pin the
    # tiering logic.
    "best_neighbour_uniformity", "neighbour_ratio",
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
          rules: Path | None = None,
          extra: list[str] | None = None) -> dict[str, dict]:
    """Run score_evidence.py and return findings keyed by sample."""
    out = tmp / f"findings_{cov.stem}.tsv"
    cmd = [sys.executable, BIN / "score_evidence.py", "--coverage", cov,
           "--watchlist", PACK / "watchlist.resolved.tsv",
           "--pack-yaml", PACK / "pack.yaml", "--out", out]
    if markers:
        cmd += ["--marker-evidence", markers,
                "--marker-rules", rules or (PACK / "marker_rules.yaml")]
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
        ("selftest:mask_pack_refs", [sys.executable, SCRIPTS / "mask_pack_refs.py", "--self-test"]),
        ("selftest:measure_reference_span",
         [sys.executable, SCRIPTS / "measure_reference_span.py", "--self-test"]),
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

    # No spaces: the SAM spec forbids them in @SQ SN, and samtools' BED reader splits on
    # whitespace, so a spaced name cannot be addressed as a region. The fixture carried
    # FASTA-style descriptions here, which made it an invalid header that happened not to
    # matter until coverage_stats.py began restricting `samtools depth` to the target's
    # contigs. Real bundle references are `taxid|<taxid>|<accession>`.
    tgt = "taxid|632|NZ_TEST1.1"
    nb = "taxid|633|NC_TEST2.1"
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
    p = run([sys.executable, BIN / "resolve_watchlist.py",
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
        # Not passed --expected-arms, so no absent arm may be called a failure.
        record("screen:an arm nobody asked for is `not_run`, never a failure",
               all(s.get(a, {}).get("status") == "not_run"
                   and not s.get(a, {}).get("flag")
                   for a in ("metabuli", "sylph", "bracken")),
               str({a: s.get(a, {}).get("status") for a in
                    ("metabuli", "sylph", "bracken")}))
    else:
        record("screen:unclassified fraction is reported and flagged", False, "no summary")
        record("screen:an arm nobody asked for is `not_run`, never a failure",
               False, "no summary")

    # AN ENABLED ARM THAT PRODUCED NOTHING IS A HOLE IN COVERAGE, NOT A NEGATIVE.
    #
    # METABULI, KRAKENUNIQ and SYLPH run under errorStrategy 'ignore' so one dead arm
    # cannot sink the run, and fillMissing() then pads the crash to the same `[]` a
    # DISABLED arm gets. By the time parse_screen.py runs, "skipped", "crashed" and
    # "emitted an empty report" are one state: no argument. Only the workflow knows
    # which it was, so it passes --expected-arms - and if that is ever dropped from
    # the process the difference silently disappears again, which is what this pins.
    summ2 = tmp / "summary_dead.tsv"
    empty_mb = tmp / "mb_empty.tsv"
    empty_mb.write_text("")
    p = run([sys.executable, BIN / "parse_screen.py", "--kraken2-report", k2,
             "--metabuli-report", empty_mb,
             "--expected-arms", "kraken2,bracken,metabuli,sylph",
             "--watchlist", wl, "--descendants", desc,
             "--sample", "S2", "--out", tmp / "cand2.tsv", "--out-summary", summ2])
    if summ2.exists():
        with summ2.open() as fh:
            s2 = {r["classifier"]: r for r in csv.DictReader(fh, delimiter="\t")}
        record("screen:an enabled arm that produced nothing is named, not absent",
               s2.get("sylph", {}).get("status") == "no_output"
               and s2.get("sylph", {}).get("flag") == "ARM_PRODUCED_NOTHING",
               str(s2.get("sylph", {})))
        record("screen:an enabled arm's empty report is not a clean negative",
               s2.get("metabuli", {}).get("status") == "empty_report"
               and s2.get("metabuli", {}).get("flag") == "ARM_EMPTY_REPORT",
               str(s2.get("metabuli", {})))
        # bracken cannot raise a candidate and exits non-zero by design on a report
        # with no reads at the requested level - which is a CLEAN NEGATIVE CONTROL.
        # Flagging that would put a warning on every spotless control in the run.
        record("screen:bracken missing is recorded but not flagged",
               s2.get("bracken", {}).get("status") == "no_output"
               and not s2.get("bracken", {}).get("flag"),
               str(s2.get("bracken", {})))
        # `0` unclassified out of `0` total is a claim about the data; sylph has no
        # read ledger at all and bracken re-estimates over kraken2's, so reporting
        # either would read as corroboration from an independent arm.
        record("screen:an arm with no read ledger reports blank, not zero",
               all(s2.get(a, {}).get("unclassified_fraction") == ""
                   and s2.get(a, {}).get("total_reads") == ""
                   for a in ("sylph", "bracken")),
               str({a: s2.get(a, {}) for a in ("sylph", "bracken")}))
        record("screen:every arm gets a row, so coverage is readable as a table",
               set(s2) == {"krakenuniq", "kraken2", "bracken", "metabuli", "sylph"},
               str(sorted(s2)))
    else:
        for t in ("an enabled arm that produced nothing is named, not absent",
                  "an enabled arm's empty report is not a clean negative",
                  "bracken missing is recorded but not flagged",
                  "an arm with no read ledger reports blank, not zero",
                  "every arm gets a row, so coverage is readable as a table"):
            record(f"screen:{t}", False, "no summary written")
    # AN EMPTY REPORT IS ONLY AN ALARM IF THE ARM SEARCHED EVERYTHING.
    #
    # Found on the first real Illumina re-score: sylph profiles against the PACK -
    # 371 watchlist references, nothing else - and emits one row per containment hit,
    # so zero rows is how it says "no watchlist organism here". Four of eleven
    # single-organism libraries were isolates of organisms nobody is watching for, and
    # all four were flagged ARM_EMPTY_REPORT and bannered at the top of the report as
    # a classifier that had produced nothing. The pipeline had in fact answered
    # correctly. kraken2, krakenuniq and metabuli are the opposite case: parsed
    # unfiltered off a comprehensive database, they emit rows for any input at all, so
    # zero usable rows there really is a broken report. The two must not share a flag,
    # because a warning that fires on correct behaviour is one the reader learns to
    # skip - and ARM_PRODUCED_NOTHING, which is never routine, is in the same column.
    summ3 = tmp / "summary_packarm.tsv"
    empty_sy = tmp / "sy_empty.tsv"
    empty_sy.write_text("Sample_file\tGenome_file\tTaxonomic_abundance\tAdjusted_ANI\n")
    run([sys.executable, BIN / "parse_screen.py", "--kraken2-report", k2,
         "--sylph-profile", empty_sy,
         "--expected-arms", "kraken2,sylph",
         "--watchlist", wl, "--descendants", desc,
         "--sample", "S3", "--out", tmp / "cand3.tsv", "--out-summary", summ3])
    if summ3.exists():
        with summ3.open() as fh:
            s3 = {r["classifier"]: r for r in csv.DictReader(fh, delimiter="\t")}
        record("screen:a pack-restricted arm's empty report is the negative answer",
               s3.get("sylph", {}).get("status") == "empty_report"
               and not s3.get("sylph", {}).get("flag"),
               str(s3.get("sylph", {})))
        record("screen:and it still says so, rather than going quiet",
               "negative answer" in s3.get("sylph", {}).get("note", ""),
               "suppressing the flag must not suppress the explanation; the row is "
               "still what tells a reader the arm ran at all")
        record("screen:a crashed pack-restricted arm is STILL flagged",
               s3.get("sylph", {}).get("status") == "empty_report"
               and _crashed_sylph_flag(tmp, k2, wl, desc) == "ARM_PRODUCED_NOTHING",
               "an absent file means the arm died; only a present-but-empty one is "
               "the negative answer, and the suppression must not reach the crash")
    else:
        for name in ("screen:a pack-restricted arm's empty report is the negative answer",
                     "screen:and it still says so, rather than going quiet",
                     "screen:a crashed pack-restricted arm is STILL flagged"):
            record(name, False, "no summary")

    # A summary is written even when EVERY arm died, which is the run it has the most
    # to say about. The writer used to be guarded by `if summary`, so that run
    # produced no file at all and the gap was indistinguishable from a missing step.
    summ3 = tmp / "summary_allfail.tsv"
    run([sys.executable, BIN / "parse_screen.py",
         "--expected-arms", "kraken2,metabuli,sylph",
         "--watchlist", wl, "--descendants", desc,
         "--sample", "S3", "--out", tmp / "cand3.tsv", "--out-summary", summ3])
    ok = False
    if summ3.exists():
        with summ3.open() as fh:
            s3 = [r for r in csv.DictReader(fh, delimiter="\t")]
        ok = (len(s3) == 5
              and {r["classifier"] for r in s3 if r["flag"] == "ARM_PRODUCED_NOTHING"}
              == {"kraken2", "metabuli", "sylph"})
    record("screen:a run where every arm died still writes its coverage table", ok,
           "no file" if not summ3.exists() else "")

    # `not_applicable` IS NOT A SHADE OF `not_run`. Found on the first real ONT run,
    # where bracken - which cannot run on long reads, its model taking one fixed read
    # length - was reported "disabled for this run" on a sample whose kraken2 row read
    # `ok` and where nobody had disabled anything. A reader would conclude an operator
    # passed --skip_kraken2. An operator's choice and a property of the platform are
    # different answers to why something was not searched.
    summ5 = tmp / "summary_ont.tsv"
    run([sys.executable, BIN / "parse_screen.py", "--kraken2-report", k2,
         "--expected-arms", "kraken2", "--inapplicable-arms", "bracken",
         "--platform", "ont", "--watchlist", wl, "--descendants", desc,
         "--sample", "S5", "--out", tmp / "cand5.tsv", "--out-summary", summ5])
    ok = False
    if summ5.exists():
        with summ5.open() as fh:
            s5 = {r["classifier"]: r for r in csv.DictReader(fh, delimiter="\t")}
        ok = (s5.get("bracken", {}).get("status") == "not_applicable"
              and "platform" in s5.get("bracken", {}).get("note", "")
              and s5.get("metabuli", {}).get("status") == "not_run")
    record("screen:an arm the platform rules out is not reported as switched off", ok,
           str(s5.get("bracken", {})) if summ5.exists() else "no summary")

    # Naming an arm as both expected and inapplicable is a contradiction in the
    # workflow's own description of the run, so it must not resolve to one of them.
    p = run([sys.executable, BIN / "parse_screen.py", "--kraken2-report", k2,
             "--expected-arms", "kraken2,bracken", "--inapplicable-arms", "bracken",
             "--watchlist", wl, "--descendants", desc, "--sample", "S6",
             "--out", tmp / "cand6.tsv", "--out-summary", tmp / "summary6.tsv"])
    record("screen:an arm cannot be both expected and inapplicable",
           p.returncode != 0 and "bracken" in p.stderr, p.stderr.strip()[-120:])

    # An unknown arm name must fail loudly: it would otherwise be a silently
    # unenforced expectation, which is the same shape as a misspelled assertion.
    p = run([sys.executable, BIN / "parse_screen.py", "--kraken2-report", k2,
             "--expected-arms", "kraken2,krakn2", "--watchlist", wl,
             "--descendants", desc, "--sample", "S4", "--out", tmp / "cand4.tsv",
             "--out-summary", tmp / "summary4.tsv"])
    record("screen:a misspelled expected arm is rejected, not ignored",
           p.returncode != 0 and "krakn2" in p.stderr, p.stderr.strip()[-120:])


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

        # THE SAME ORPHAN CHECK FOR THE SECOND EXPECTATION AXIS.
        #
        # A `- marker:` expectation naming an id that is not in markers.tsv is an
        # expectation that can never pass - and the lesson this file already
        # records for check_tier_name is that a misspelled assertion silently
        # removes the gate it was written to be. The evaluator reports it as
        # ERROR at run time; catching it statically means a typo fails the suite
        # instead of a benchmark run hours later.
        with (PACK / "markers.tsv").open() as fh:
            marker_ids = {r["marker_id"] for r in csv.DictReader(
                [l for l in fh if not l.startswith("#")], delimiter="\t")
                if r.get("marker_id")}
        want_markers, bad_states = set(), []
        for sn in sc:
            for e in sn.get("expect") or []:
                if "marker" not in e:
                    continue
                want_markers.add(e["marker"])
                if str(e.get("state", "")).strip().upper() not in (
                        "PRESENT", "ABSENT", "UNKNOWN"):
                    bad_states.append(f"{sn.get('id')}/{e['marker']}={e.get('state')!r}")
        record("pack:every scenario marker expectation names a real marker",
               not (want_markers - marker_ids),
               f"orphans: {sorted(want_markers - marker_ids)}")
        record("pack:every scenario marker expectation names a real state",
               not bad_states, f"bad: {bad_states}")

        # EVERY CATALOGUE ROW MUST BE ABLE TO BUILD A QUERY.
        #
        # query_names() drops names under MIN_QUERY_NAME_LEN, so a row whose gene
        # symbol is one character with no alias produces NO query and is never
        # searched - henipa_cedar (gene 'N', aliases '-') sat that way, recorded as
        # `failed`, which reads as "NCBI has nothing" rather than "we never asked".
        # Its veto for BOTH Nipah and Hendra could therefore never fire. Eighteen
        # sibling rows carry a one-letter gene symbol rescued by an alias; this is
        # the check that a nineteenth cannot be added without one.
        try:
            sys.path.insert(0, str(SCRIPTS))
            import build_marker_db as bmd
            unsearchable = [
                r["marker_id"] for r in bmd.read_marker_catalogue(PACK / "markers.tsv")
                if r["query_strategy"] not in ("curate_required", "motif")
                and not bmd.nt_queries(r) and not bmd.aa_queries(r)
            ]
            record("pack:every retrievable marker can build at least one query",
                   not unsearchable, f"build no query: {unsearchable}")
        except Exception as exc:  # noqa: BLE001
            record("pack:every retrievable marker can build at least one query",
                   False, f"{type(exc).__name__}: {exc}")
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

        # ---- read/breadth gate coherence, PINNED ----
        #
        # The two stage-2 gates are ANDed, and breadth is bounded by read count
        # (best case N*L/G), so whenever min_breadth*G/L > min_unique_reads the read
        # floor cannot bind and the watchlist's stated LoD understates the real one by
        # that ratio - 243x for Schistosoma mansoni as shipped. These numbers are due
        # to be set from the benchmark, so rather than recalibrate on a made-up spread
        # factor, the counts are pinned here: a pack edit then shows up as coherence
        # improving or regressing instead of silently staying broken.
        import subprocess
        refs_manifest = Path(os.environ.get(
            "BTF_REFS_MANIFEST",
            "/media/8tb/bioThreatFinder/packs/human/refs/manifest.tsv"))
        if not refs_manifest.is_file():
            skip("pack:read/breadth gate coherence is tracked",
                 "no built refs manifest (set BTF_REFS_MANIFEST)")
        else:
            cp = subprocess.run(
                [sys.executable, str(SCRIPTS / "check_threshold_coherence.py"),
                 "--refs-manifest", str(refs_manifest)],
                capture_output=True, text=True)
            out = cp.stdout
            def _num(label: str) -> int:
                for ln in out.splitlines():
                    if ln.startswith(label):
                        return int(ln.split(":", 1)[1].strip().split()[0])
                return -1
            assessed, breadth_binds = _num("taxa assessed"), _num("breadth is the real gate")
            # KNOWN DEBT, recorded so it cannot drift unnoticed. Lowering this number is
            # an improvement; raising it means a pack edit made the mismatch worse.
            # 89 -> 90 and 75 -> 76 when the Filoviridae group entry was added (D25).
            # It inherits the filovirus species' 10-read floor against a 220 kb
            # 12-genome bundle, so breadth binds first, exactly as for the species
            # rows it sits above. One more instance of known debt, not a new kind.
            record("pack:read/breadth gate coherence is tracked",
                   assessed == 90 and breadth_binds == 76,
                   f"{breadth_binds} of {assessed} taxa have a decorative read floor "
                   f"(baseline 75 of 89). See docs/decisions.md D11.")
    except Exception as exc:
        record("pack:scenario file parses", False, f"{type(exc).__name__}: {exc}")


# ---------------------------------------------------------------------------
# 7b. Depth-conditional strain / sub-species typing
# ---------------------------------------------------------------------------

def test_depth_gate(tmp: Path) -> None:
    """A strain claim must be refused when the depth cannot support it.

    Marker typing survives low depth (short locus, not a whole genome) but ABSENCE does
    not: a single-copy 1.5 kb marker at 0.012x expects 0.12 reads, so finding none is
    the expected outcome whether or not the locus is there. Reporting ABSENT there turns
    no data into positive evidence against a real finding.
    """
    sys.path.insert(0, str(BIN))
    import score_evidence as se

    th = {"min_expected_reads_for_absence": 3.0, "assumed_marker_bp": 1000.0,
          "min_depth_for_genome_typing": 20.0}

    # The arithmetic itself: reads = depth * length / read_len.
    record("depth:expected reads at 0.62x on a 1.5 kb marker is ~6",
           abs(se.expected_marker_reads(0.62, 1500, 150) - 6.2) < 0.01,
           f"{se.expected_marker_reads(0.62, 1500, 150):.2f}")
    record("depth:expected reads at 0.012x on a 1.5 kb marker is ~0.12",
           se.expected_marker_reads(0.012, 1500, 150) < 0.2,
           f"{se.expected_marker_reads(0.012, 1500, 150):.3f}")

    try:
        import marker_rules as mr
    except Exception as exc:
        record("depth:marker_rules importable", False, str(exc))
        return

    metas = {
        "ecoli_stx1": {"taxon": "Escherichia coli", "marker_class": "pathotype",
                       "nt_length": 1500.0, "nt_reads": 0.0, "state": "ABSENT"},
        "other_taxon": {"taxon": "Yersinia pestis", "marker_class": "pathotype",
                        "nt_length": 1500.0, "nt_reads": 0.0, "state": "ABSENT"},
    }
    states = {"ecoli_stx1": mr.FALSE, "other_taxon": mr.FALSE}

    # Low depth: the absence is not interpretable, so it must become UNKNOWN.
    adj, demoted = se.gate_absence_on_depth(states, metas, "Escherichia coli",
                                            0.012, 150.0, th)
    record("depth:an uninterpretable absence becomes UNKNOWN",
           "ecoli_stx1" not in adj and len(demoted) == 1,
           f"adjusted={adj}, demoted={demoted}")
    record("depth:another taxon's marker is left alone",
           adj.get("other_taxon") == mr.FALSE,
           "depth is per taxon; a neighbour's marker must not be rejudged")

    # Ample depth: the absence IS evidence and must survive.
    adj2, demoted2 = se.gate_absence_on_depth(states, metas, "Escherichia coli",
                                              30.0, 150.0, th)
    record("depth:a real absence at good depth stays ABSENT",
           adj2.get("ecoli_stx1") == mr.FALSE and not demoted2,
           f"adjusted={adj2}, demoted={demoted2}")

    # ---- depth-aware marker PRESENCE, the mirror of the above ----
    #
    # A fixed min_coverage of 0.80 silently demands >= 1.6x on the locus. Below that it
    # is unreachable however real the organism is. Measured on the plumbing control:
    # ypest_pmt1_caf1, 13 reads at 99.8% identity, 0.323 breadth against a 0.80 floor,
    # which held a correct Y. pestis detection at screened_only.
    import marker_evidence as me
    record("depth:achievable breadth at 0.34x is ~0.29, so a 0.80 floor is unreachable",
           abs(me.expected_breadth_at(0.34) - 0.288) < 0.005
           and me.expected_breadth_at(0.34) < 0.80,
           f"{me.expected_breadth_at(0.34):.3f}")
    record("depth:a 0.80 floor becomes reachable by ~3x",
           me.expected_breadth_at(3.0) > 0.80, f"{me.expected_breadth_at(3.0):.3f}")
    record("depth:expected breadth is 0 at zero depth",
           me.expected_breadth_at(0.0) == 0.0)
    # The caf1 case: observed breadth clears 0.60 x achievable.
    caf1_depth, caf1_obs = 0.47, 0.323
    record("depth:caf1 qualifies on the depth-limited route",
           caf1_obs >= 0.60 * me.expected_breadth_at(caf1_depth),
           f"observed {caf1_obs} vs required "
           f"{0.60 * me.expected_breadth_at(caf1_depth):.3f}")
    # And it must NOT qualify on the strict route, or the relaxation is meaningless.
    record("depth:caf1 still fails the strict 0.80 floor",
           caf1_obs < 0.80, "the depth-limited route is what admits it, not the floor")

    # And the genome-wide claim must say it was not attempted.
    claim, why = se.describe_strain_claim(0.62, th)
    record("depth:genome typing is reported NOT ATTEMPTED at 0.62x",
           claim == "marker_typing_only" and "NOT ATTEMPTED" in why, f"{claim}: {why[:70]}")
    claim2, why2 = se.describe_strain_claim(25.0, th)
    record("depth:genome typing is flagged supportable at 25x",
           claim2 == "genome_typing_supported", f"{claim2}")
    record("depth:and says it is not implemented yet, rather than implying a result",
           "NOT YET IMPLEMENTED" in why2, why2[:80])


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


def test_workflow_scripts(tmp: Path) -> None:
    """Every script the WORKFLOW invokes must be in bin/ and executable.

    Nextflow puts only `bin/` on the task PATH. `resolve_watchlist.py` lived in
    `scripts/`, so RESOLVE_WATCHLIST died with exit 127 (command not found) on the very
    first end-to-end run - which means the pipeline had never been runnable at all. It
    failed loudly, which is the good case, but it cost a full launch to discover and
    nothing in the suite would have caught it.
    """
    import re
    nf_text = "\n".join(
        f.read_text(encoding="utf-8")
        for f in [REPO / "main.nf", *sorted((REPO / "modules").rglob("*.nf"))])

    # Script invocations: a .py at the start of a command line inside a script block.
    invoked = sorted(set(re.findall(r"^\s{4,}([a-z_][a-z0-9_]*\.py)\s", nf_text,
                                    re.MULTILINE)))
    record("workflow:at least one script invocation was found to check",
           len(invoked) >= 5, f"found {invoked}")

    missing, not_exec = [], []
    for name in invoked:
        target = BIN / name
        if not target.is_file():
            missing.append(name)
        elif not os.access(target, os.X_OK):
            not_exec.append(name)
    record("workflow:every invoked script is in bin/ (Nextflow's task PATH)",
           not missing,
           f"not in bin/: {missing} - these would fail with exit 127")
    record("workflow:every invoked script is executable",
           not not_exec, f"not executable: {not_exec}")

    # ---- optional inputs must not use a shared placeholder file ----
    #
    # fillMissing() used to pad every absent stage-1 arm with one shared assets/NO_FILE.
    # Nextflow stages inputs by basename, so a process receiving two or more of them died
    # with "input file name collision" - which PARSE_SCREEN did the moment two arms were
    # skipped. An absent optional path is now an empty list: nothing is staged, so there
    # is no basename that can collide. This asserts the placeholder has not crept back,
    # which is a stronger and simpler guarantee than checking for stageAs workarounds.
    nf_files = [REPO / "main.nf", *sorted((REPO / "modules").rglob("*.nf"))]
    offenders = []
    for nf_file in nf_files:
        for i, line in enumerate(nf_file.read_text(encoding="utf-8").splitlines(), 1):
            code = re.sub(r"//.*$", "", line)          # ignore comments
            if "NO_FILE" in code:
                offenders.append(f"{nf_file.name}:{i}")
    record("workflow:no shared NO_FILE placeholder is used for optional inputs",
           not offenders,
           f"found at {offenders[:6]} - an absent optional path should be [], which "
           f"stages nothing and cannot collide")
    record("workflow:the NO_FILE asset is gone",
           not (REPO / "assets" / "NO_FILE").exists(),
           "assets/NO_FILE still present")

    # ---- stage 2 must not re-gate candidates on a read count ----
    #
    # main.nf carried a
    #     .filter { (row.assigned_reads as Integer) >= params.min_screen_reads }
    # between PARSE_SCREEN and stage 2. `assigned_reads` is a READ COUNT and only the
    # read-based arms produce one; sylph reports ANI containment, so a sylph-raised
    # candidate has assigned_reads=0 by construction and was silently dropped. On the
    # first end-to-end run that discarded Y. pestis (sylph abundance 6.46) and E. coli
    # (22.17) and the run exited 0 with a header-only findings.tsv - a clean negative
    # for the whole watchlist, reported as success, on a positive control.
    #
    # bin/parse_screen.py is the single place candidacy is decided and it already
    # applies --min-reads per classifier.
    # Strip // comments first: the comment explaining this defect quotes the removed
    # filter verbatim, and matching that would make the test fail on its own docs.
    main_text = "\n".join(
        re.sub(r"//.*$", "", ln)
        for ln in (REPO / "main.nf").read_text(encoding="utf-8").splitlines())
    gate = re.search(r"\.filter\s*\{[^}]*assigned_reads[^}]*\}", main_text)
    record("workflow:stage 2 does not re-filter candidates on assigned_reads",
           gate is None,
           f"found a read-count gate after PARSE_SCREEN: "
           f"{gate.group(0)[:90] if gate else ''} - this drops every candidate raised "
           f"by a classifier that reports containment rather than reads")


def test_shared_assembly_masking(tmp: Path) -> None:
    """Taxa sharing one NCBI assembly must carry IDENTICAL masking numbers.

    F. tularensis (263) / F. novicida (264) share GCF_000833355.1, and
    S. dysenteriae (622) / Shigella (620) share GCF_022354085.1. Because the FASTAs
    carry identical contig names, mask_pack_refs.py used to attribute the contigs to
    whichever genome loaded first and skip the second, writing genome_bp=0 masked_bp=0
    and an interval-free BED - byte-identical in meaning to "examined and found clean".
    Both shipped BEDs were in that state until 2026-10-03.

    The invariant is agreement: identical sequence must produce identical masking.
    The SHARED_SEQUENCE flag is now asserted too - it used to be set from RUN-TIME
    contig collisions, so it only appeared when both partners were in the same
    invocation and a --only-taxid run of one partner produced correct numbers with an
    empty flag. It is now derived from the refs manifest's `accession` column, which
    is an artifact of the pack rather than of how the script was called.

    Skips rather than fails when the built artifacts are absent, since the suite is
    meant to run offline on a checkout without /media/8tb.
    """
    import csv as _csv
    refs_manifest = pathlib.Path(
        "/media/8tb/bioThreatFinder/packs/human/refs/manifest.tsv")
    mask_manifest = pathlib.Path(
        "/media/8tb/bioThreatFinder/packs/human/refs/masked/mask_manifest.tsv")
    if not mask_manifest.exists():
        skip("shared_assembly:pairs agree", f"no built artifacts at {mask_manifest}")
        return
    # KEYED ON THE REFERENCE, NOT THE TAXID. This dict used to be {taxid: row}, and
    # taxid 11266 has TWELVE rows, so eleven of them were silently discarded by the
    # comprehension - the same taxid-as-key mistake that made the manifest writer
    # drop rows it was not rewriting. The pairs checked below happen to be
    # one-reference taxa, so the test passed; it was asserting less than it looked.
    with mask_manifest.open() as fh:
        mask_rows = list(_csv.DictReader(fh, delimiter="\t"))
    rows: dict = {}
    for r in mask_rows:
        rows.setdefault(r["taxid"], r)
    record("shared_assembly:no reference is lost to taxid-keyed bookkeeping",
           len({r.get("bed_path") for r in mask_rows}) == len(mask_rows),
           f"{len(mask_rows)} rows, "
           f"{len({r.get('bed_path') for r in mask_rows})} distinct bed_path")

    for a, b, label in (("263", "264", "F. tularensis / F. novicida"),
                        ("620", "622", "Shigella / S. dysenteriae")):
        ra, rb = rows.get(a), rows.get(b)
        if not ra or not rb:
            skip(f"shared_assembly:{label}", f"taxid {a if not ra else b} not in manifest")
            continue
        # The PER-SOURCE columns are the point. The first version of this test compared
        # only the four totals, all of which agreed, so it passed green while
        # crossgroup_bp read 1,835 for taxid 620 and 2,204 for 622 - two rows whose BED
        # intervals are identical. A guard whose field list omits the columns that
        # disagree asserts nothing.
        fields = ("genome_bp", "masked_bp", "masked_fraction", "n_masked_intervals",
                  "lowcomplexity_bp", "vector_bp", "host_bp", "crossgroup_bp")
        diffs = [f"{f}: {ra.get(f)} vs {rb.get(f)}"
                 for f in fields if ra.get(f) != rb.get(f)]
        record(f"shared_assembly:{label} agree", not diffs, "; ".join(diffs))
        # The vacuous state the bug produced, asserted directly so a regression that
        # reintroduces it cannot pass by agreeing at zero.
        record(f"shared_assembly:{label} not vacuous",
               ra.get("genome_bp") not in ("", "0") and rb.get("genome_bp") not in ("", "0"),
               f"genome_bp {ra.get('genome_bp')} / {rb.get('genome_bp')}")
        # Now assertable, because the flag no longer depends on how the script was
        # invoked. Each partner must NAME the other: agreement makes the numbers
        # auditable, the flag is what says the agreement is intentional.
        record(f"shared_assembly:{label} say so in the flag",
               "SHARED_SEQUENCE" in (ra.get("flag") or "")
               and "SHARED_SEQUENCE" in (rb.get("flag") or ""),
               f"flags {ra.get('flag')!r} / {rb.get('flag')!r}")

    # EVERY PACK REFERENCE MUST HAVE A MASKING RESULT.
    #
    # The invariant above is about two references agreeing. This one is about a
    # reference existing at all, and it is the one that was being violated: on
    # 2026-10-15 the pack held 371 references and 359 mask BEDs, with all twelve of
    # the Filoviridae group entry's files unmasked because the group entry (D24/D25)
    # was added AFTER the masking pass. Six of the twelve had contigs in no BED
    # anywhere, so coverage_stats scored them against a mask map that simply did not
    # contain them - indistinguishable from examined-and-clean. Found by hand-diffing
    # a directory listing against a manifest, which is why it is asserted here and
    # reconciled by the script itself.
    if refs_manifest.exists():
        with refs_manifest.open() as fh:
            refs = [r for r in _csv.DictReader(fh, delimiter="\t")
                    if (r.get("path") or "").strip()]
        have = {pathlib.Path(r["bed_path"]).name for r in mask_rows if r.get("bed_path")}
        missing = [r for r in refs
                   if pathlib.Path(r["path"]).name
                   .replace(".fna.gz", "").replace(".fna", "") + ".mask.bed" not in have]
        record("shared_assembly:every pack reference has a masking result",
               not missing,
               f"{len(missing)} unmasked, e.g. "
               + ", ".join(sorted({r['name'] for r in missing})[:3]))
    else:
        skip("shared_assembly:every pack reference has a masking result",
             "no refs manifest")


def test_exhaustive_absent(tmp: Path) -> None:
    """truth_absent must be the whole watchlist minus what is really in the sample.

    The evaluator scores a call for an unlisted taxon as `unlisted`, not
    `false_positive`, so a hand-written absent list silently excludes from the
    specificity table any organism its author forgot. negative_host_only listed eleven
    taxa, omitted Toxoplasma gondii, and the 2026-10-03 run reported T. gondii
    `supported` in all 21 samples off 100% human reads while specificity read 0.000 with
    zero T. gondii rows. docs/decisions.md D14.

    Asserts both halves: the forgotten organism IS now absent, and the things genuinely
    in the sample are NOT - because listing a present organism as absent would score a
    CORRECT detection as a false positive, which is the opposite failure.
    """
    sys.path.insert(0, str(REPO / "scripts"))
    try:
        import simulate_benchmark as sb
    except Exception as exc:                                    # noqa: BLE001
        skip("absent:exhaustive", f"cannot import simulate_benchmark: {exc}")
        return

    class _Spike:
        def __init__(self, taxon): self.taxon = taxon

    class _Comp:
        def __init__(self, label, kind): self.label, self.kind = label, kind

    class _Sc:
        spikes = [_Spike("Yersinia pestis")]
        vaccine = ["Mycobacterium bovis"]
        absent = ["Bacillus anthracis"]

    watchlist = ["Bacillus anthracis", "Yersinia pestis", "Toxoplasma gondii",
                 "Escherichia coli", "Mycobacterium bovis", "Francisella tularensis"]
    components = [_Comp("Escherichia coli", "background"),
                  _Comp("Yersinia pestis", "spike")]
    got = sb.exhaustive_absent(_Sc(), ["Yersinia pestis"], watchlist, components)
    low = {g.lower() for g in got}

    record("absent:forgotten watchlist organism is now absent",
           "toxoplasma gondii" in low, f"got {sorted(got)}")
    record("absent:unrelated watchlist organism is absent",
           "francisella tularensis" in low, f"got {sorted(got)}")
    record("absent:explicit entry is kept",
           "bacillus anthracis" in low, f"got {sorted(got)}")
    # Each of these would turn a correct detection into a scored false positive.
    record("absent:spiked target is NOT absent",
           "yersinia pestis" not in low, f"got {sorted(got)}")
    record("absent:background organism is NOT absent",
           "escherichia coli" not in low, f"got {sorted(got)}")
    record("absent:vaccine strain is NOT absent",
           "mycobacterium bovis" not in low, f"got {sorted(got)}")


def test_matrix_blank_suppresses_host_artifact(tmp: Path) -> None:
    """A host-derived false positive is only suppressible if a control CONTAINS host.

    D14: *Toxoplasma gondii* was reported `supported` in 18/18 samples of lod_v1 and 6/6
    of divergent_v1, off reads that attribution showed were 100% human. The suppression
    machinery was never broken - the batch's only control was an NTC built with
    host_fraction 0.0, so a human-derived artefact could not appear in it, control_rpm
    stayed 0 for that taxon, and `suppressed_control` could never fire.

    Both halves are asserted, because only the pair distinguishes the fix from a tier
    that happens to be right. With an NTC that has no such row the finding must remain
    `supported` - reproducing the bug - and with a matrix blank carrying the artefact at
    a comparable level it must become `suppressed_control`. A test that checked only the
    second would also pass if everything were suppressed unconditionally.
    """
    shape = dict(taxid="5811", taxon="Toxoplasma gondii", group="parasite",
                 reference_bp=65000000, breadth=0.4410, mean_depth=0.02,
                 window_bp=1000, windows_expected=65000, windows_occupied=62725,
                 window_occupancy=0.9650, mean_identity=0.9815,
                 neighbour_ratio="inf", mapped_fraction_of_assigned="1.0000",
                 total_reads_screened=2000000, mapq_min=10)

    # (a) NTC present but, having no host, carrying no T. gondii row at all.
    cov_a = tmp / "blank_absent.tsv"
    write_cov(cov_a, [
        cov_row(sample="SPECIMEN", role="sample", unique_reads=8310,
                reads_per_million=4155.0, depth_cv=5.44, screen_assigned_reads=8400,
                **shape),
        cov_row(sample="NTC", role="control", taxid="562", taxon="Escherichia coli",
                group="bacteria", reference_bp=4600000, unique_reads=900,
                reads_per_million=4500.0, breadth=0.2000, mean_depth=0.03,
                depth_cv=1.20, window_bp=1000, windows_expected=4600,
                windows_occupied=920, window_occupancy=0.2000, mean_identity=0.9990,
                neighbour_ratio="inf", screen_assigned_reads=950,
                mapped_fraction_of_assigned="0.9474",
                total_reads_screened=200000, mapq_min=10),
    ])
    got_a = score(tmp, cov_a).get("SPECIMEN", {}).get("tier", "<missing>")
    record("matrix_blank:host artefact is NOT suppressed by a host-free NTC",
           got_a == "supported", f"got {got_a}, want supported (the D14 bug)")

    # (b) Same specimen, plus a matrix blank that DOES carry the artefact.
    cov_b = tmp / "blank_present.tsv"
    write_cov(cov_b, [
        cov_row(sample="SPECIMEN", role="sample", unique_reads=8310,
                reads_per_million=4155.0, depth_cv=5.44, screen_assigned_reads=8400,
                **shape),
        cov_row(sample="BLANK", role="control", unique_reads=8447,
                reads_per_million=4267.0, depth_cv=5.50, screen_assigned_reads=8500,
                **shape),
    ])
    f_b = score(tmp, cov_b)
    got_b = f_b.get("SPECIMEN", {}).get("tier", "<missing>")
    record("matrix_blank:host artefact IS suppressed by a matrix blank",
           got_b == "suppressed_control", f"got {got_b}, want suppressed_control")
    # The blank's own content is QC, never a finding about a specimen.
    got_c = f_b.get("BLANK", {}).get("tier", "<missing>")
    record("matrix_blank:the blank's own row is control_background",
           got_c == "control_background", f"got {got_c}")


def test_platform_carried(tmp: Path) -> None:
    """Every finding must carry ITS OWN row's platform, whatever tier it landed in.

    `plat` was read only inside the innermost else, so a row that short-circuited
    earlier (read floor, neighbour attribution, control suppression) fell back on
    whatever Python had left in the loop local. First row short-circuiting =>
    UnboundLocalError and no findings at all; any later row => the PREVIOUS row's
    platform silently written into this finding. That field tells a reader which
    read-identity floor was applied, so a stale value misreports an ONT confirmation
    as Illumina. Found by the 2026-10-03 Kraken2 LoD run, which crashed here.

    The ordering matters: the FIRST row must be one that short-circuits.
    """
    cov = tmp / "platform_carried.tsv"
    write_cov(cov, [
        # Below its read floor -> stops before the identity-floor branch. First on purpose.
        cov_row(sample="ONT_LOWREADS", role="sample", taxid="632",
                taxon="Yersinia pestis", group="bacteria", platform="ont",
                reference_bp=4650000, unique_reads=1, reads_per_million=0.2,
                breadth=0.0001, mean_depth=0.01, depth_cv=9.0, window_bp=1000,
                windows_expected=4650, windows_occupied=1, window_occupancy=0.0002,
                mean_identity=0.9100, total_reads_screened=5000000, mapq_min=10),
        # Reaches the full gate set, so it assigns plat itself.
        cov_row(sample="ILL_FULL", role="sample", taxid="632", taxon="Yersinia pestis",
                group="bacteria", platform="illumina", reference_bp=4650000,
                unique_reads=820, reads_per_million=164.0, breadth=0.3400,
                mean_depth=0.88, depth_cv=1.15, window_bp=1000, windows_expected=4650,
                windows_occupied=2790, window_occupancy=0.6000, mean_identity=0.9950,
                best_neighbour_taxid="633", best_neighbour_reads=12,
                neighbour_ratio="68.33", screen_assigned_reads=1000,
                mapped_fraction_of_assigned="0.8200", total_reads_screened=5000000,
                mapq_min=10),
        # Short-circuits AFTER an ont row: the leak case, which never crashed.
        cov_row(sample="ILL_LOWREADS", role="sample", taxid="632",
                taxon="Yersinia pestis", group="bacteria", platform="illumina",
                reference_bp=4650000, unique_reads=2, reads_per_million=0.4,
                breadth=0.0002, mean_depth=0.01, depth_cv=8.0, window_bp=1000,
                windows_expected=4650, windows_occupied=2, window_occupancy=0.0004,
                mean_identity=0.9950, total_reads_screened=5000000, mapq_min=10),
    ])
    # Reaching this at all is half the test: the bug aborted before writing findings.
    got = score(tmp, cov)
    for sample, want in (("ONT_LOWREADS", "ont"), ("ILL_FULL", "illumina"),
                         ("ILL_LOWREADS", "illumina")):
        row = got.get((sample, "Yersinia pestis"))
        record(f"platform:{sample} keeps its own platform={want}",
               row is not None and row["platform"] == want,
               f"got {row and row['platform']!r}")

# ---------------------------------------------------------------------------
# 15. The attributed-to neighbour is NAMED, not just numbered
# ---------------------------------------------------------------------------

def test_neighbour_named(tmp: Path) -> None:
    """`attributed_to_neighbour` must say WHICH organism took the reads.

    This is the only tier that reports a watchlist organism as absent on the
    strength of a different organism, and it used to identify that organism by
    number alone: "reads favour neighbour taxid 620". D18 records the two calls
    that exposed it - *E. coli* attributed to 620 and *C. neoformans* to 40410,
    both organisms present at 100% abundance - where the analyst could not tell
    from the finding that 620 is *Shigella*, let alone that the pack's own notes
    say the pair cannot be separated at species level.

    The fallback is asserted as well, because the name comes from a file the
    scorer must not depend on: no manifest means the bare taxid comes back, and
    scoring itself is unchanged.
    """
    cov = tmp / "named.tsv"
    write_cov(cov, [
        cov_row(sample="NB", role="sample", taxid="632", taxon="Yersinia pestis",
                group="bacteria", reference_bp=4650000, unique_reads=300,
                reads_per_million=60.0, breadth=0.0100, mean_depth=0.10, depth_cv=9.0,
                window_bp=1000, windows_expected=4650, windows_occupied=47,
                window_occupancy=0.0100, mean_identity=0.9910,
                best_neighbour_taxid="633", best_neighbour_reads=900,
                neighbour_ratio="0.3333", best_neighbour_uniformity="0.9100",
                screen_assigned_reads=1200, mapped_fraction_of_assigned="0.7000",
                total_reads_screened=5000000, mapq_min=10),
    ])
    manifest = tmp / "manifest.tsv"
    manifest.write_text("taxid\tname\tgroup\trole\n"
                        "632\tYersinia pestis\tbacteria\ttarget\n"
                        "633\tYersinia pseudotuberculosis\tbacteria\tnear_neighbour\n")

    with_m = score(tmp / "with_m", cov, extra=["--refs-manifest", manifest])["NB"]
    record("neighbour_named:tier is attributed_to_neighbour",
           with_m["tier"] == "attributed_to_neighbour", with_m["tier"])
    record("neighbour_named:reason carries the neighbour's name",
           "Yersinia pseudotuberculosis" in with_m["reasons"],
           with_m["reasons"][:120])
    record("neighbour_named:reason still carries the taxid",
           "taxid 633" in with_m["reasons"], with_m["reasons"][:120])
    record("neighbour_named:near_neighbour role is stated",
           "near-neighbour reference" in with_m["reasons"], with_m["reasons"][:120])
    record("neighbour_named:name travels as its own column",
           with_m.get("best_neighbour_name") == "Yersinia pseudotuberculosis",
           repr(with_m.get("best_neighbour_name")))

    # No manifest: degrade to the number, and change NOTHING else.
    without = score(tmp / "without_m", cov)["NB"]
    record("neighbour_named:without a manifest the taxid is still named",
           "taxid 633" in without["reasons"], without["reasons"][:120])
    record("neighbour_named:without a manifest the name column is empty",
           without.get("best_neighbour_name") == "",
           repr(without.get("best_neighbour_name")))
    record("neighbour_named:naming does not change the tier",
           without["tier"] == with_m["tier"],
           f'{without["tier"]} vs {with_m["tier"]}')

    # The OTHER branch that names a neighbour: attribution withheld because the
    # neighbour's reads pile on shared sequence (D12). Its prose was numbered too.
    ped = tmp / "pedestal.tsv"
    write_cov(ped, [
        cov_row(sample="PED", role="sample", taxid="632", taxon="Yersinia pestis",
                group="bacteria", reference_bp=4650000, unique_reads=300,
                reads_per_million=60.0, breadth=0.0100, mean_depth=0.10, depth_cv=9.0,
                window_bp=1000, windows_expected=4650, windows_occupied=47,
                window_occupancy=0.0100, mean_identity=0.9910,
                best_neighbour_taxid="633", best_neighbour_reads=900,
                neighbour_ratio="0.3333", best_neighbour_uniformity="0.1200",
                screen_assigned_reads=1200, mapped_fraction_of_assigned="0.7000",
                total_reads_screened=5000000, mapq_min=10),
    ])
    w = score(tmp / "ped", ped,
              extra=["--refs-manifest", manifest,
                     "--min-neighbour-uniformity", "0.50"])["PED"]
    record("neighbour_named:withheld attribution is not the neighbour tier",
           w["tier"] != "attributed_to_neighbour", w["tier"])
    record("neighbour_named:withheld attribution names the neighbour",
           "Yersinia pseudotuberculosis" in w["reasons"] and "NOT CREDIBLE" in w["reasons"],
           w["reasons"][:160])


# ---------------------------------------------------------------------------
# 16. Contesting an attribution: the four cases, and the guard that matters most
# ---------------------------------------------------------------------------

def test_contest_attribution(tmp: Path) -> None:
    """A neighbour may only take a finding the target has no claim on.

    The attribution branch never looked at the target's own evidence, so a target
    covering 0.667 of its reference lost to a sibling at ratio 0.9032 and was
    reported absent - and because that sibling is not a watchlist taxon, NO finding
    anywhere carried its 3.3M reads. `--contest-attribution` asks the prior question
    first, and then the neighbour's own coverage uniformity decides which of two very
    different situations this is.

    THE NOT-CONTESTED CASE IS THE ONE THAT KEEPS THIS HONEST. D18 measured what
    happens when a neighbour gate is applied to every attributed row: 25 false
    positives on real reads. The precondition here is D14's coverage-geometry gates,
    which that measurement showed to be the discriminator (TP median breadth 0.8661
    against FP 0.0015), so a weak target's attribution is left exactly as it was -
    39 of 45 rows in zymo_rl100, 84 of 84 in lod_v4.
    """
    manifest = tmp / "cm.tsv"
    manifest.write_text("taxid\tname\tgroup\trole\n"
                        "632\tYersinia pestis\tbacteria\ttarget\n"
                        "633\tYersinia pseudotuberculosis\tbacteria\tnear_neighbour\n")

    def row(sample, role, breadth, occ, nb_unif):
        return cov_row(sample=sample, role=role, taxid="632",
                       taxon="Yersinia pestis", group="bacteria",
                       reference_bp=4650000, unique_reads=900,
                       reads_per_million=180.0, breadth=breadth, mean_depth=0.9,
                       depth_cv=1.2, window_bp=1000, windows_expected=4650,
                       windows_occupied=int(4650 * occ), window_occupancy=occ,
                       mean_identity=0.9950, best_neighbour_taxid="633",
                       best_neighbour_reads=1000, neighbour_ratio="0.9000",
                       best_neighbour_uniformity=nb_unif,
                       screen_assigned_reads=1100,
                       mapped_fraction_of_assigned="0.8200",
                       total_reads_screened=5000000, mapq_min=10)

    cov = tmp / "contest.tsv"
    write_cov(cov, [
        # Target covers its reference; neighbour's reads pile on shared sequence.
        row("PEDESTAL", "sample", 0.3400, 0.6000, "0.1500"),
        # Target covers its reference; so does the neighbour. Sibling pair.
        row("SIBLING", "sample", 0.3400, 0.6000, "0.8600"),
        # Target does NOT cover its reference. Nothing to contest.
        row("WEAK", "sample", 0.0010, 0.0200, "0.8600"),
        # Target covers its reference; the neighbour was never measured.
        row("UNMEASURED", "sample", 0.3400, 0.6000, ""),
    ])
    # --min-neighbour-uniformity 0 IS THE PIPELINE'S SETTING (D18 retracted the 0.50
    # default in nextflow.config), and it has to be passed explicitly here or the
    # comparison measures the wrong thing: at the script's own 0.50 default the
    # PEDESTAL row is already withheld by the blunt gate, so "off" would not be the
    # behaviour the pipeline actually has and the flag's effect would be understated.
    base = ["--refs-manifest", manifest, "--min-neighbour-uniformity", "0"]
    on = score(tmp / "con_on", cov,
               extra=base + ["--contest-attribution",
                             "--pedestal-max-uniformity", "0.30"])
    off = score(tmp / "con_off", cov, extra=base)

    # The control goes in its OWN batch. Sharing a batch put a 180 RPM control row on
    # taxid 632 into run-global control_rpm, so every sample row above was tested
    # against 5x its own abundance and came back suppressed_control - the fixture
    # would have asserted control subtraction while claiming to test attribution.
    ctrl = tmp / "contest_ctrl.tsv"
    write_cov(ctrl, [row("NTC", "control", 0.3400, 0.6000, "0.8600")])
    on_ctrl = score(tmp / "con_ctrl", ctrl,
                    extra=base + ["--contest-attribution",
                                  "--pedestal-max-uniformity", "0.30"])

    record("contest:pedestal neighbour cannot take the finding",
           on["PEDESTAL"]["tier"] != "attributed_to_neighbour",
           on["PEDESTAL"]["tier"])
    record("contest:pedestal rejection says why",
           "ATTRIBUTION CONTESTED AND REJECTED" in on["PEDESTAL"]["reasons"],
           on["PEDESTAL"]["reasons"][:120])
    record("contest:two covered references are indistinguishable",
           on["SIBLING"]["tier"] == "indistinguishable_complex",
           on["SIBLING"]["tier"])
    record("contest:indistinguishable names both the pair and the limit",
           "CANNOT SEPARATE" in on["SIBLING"]["reasons"]
           and "Yersinia pseudotuberculosis" in on["SIBLING"]["reasons"],
           on["SIBLING"]["reasons"][:140])
    # THE GUARD. A target that does not cover its reference has no claim to contest.
    record("contest:a weak target does NOT contest its attribution",
           on["WEAK"]["tier"] == "attributed_to_neighbour", on["WEAK"]["tier"])
    # ABSENT EVIDENCE LEAVES THE DECISION EXACTLY AS IT WAS. Not "assume pedestal"
    # (which would credit the target) and not "assume sibling" (which the first version
    # did, and which made byte-identical rows score differently in runs predating the
    # column). No measurement means nothing to contest with.
    record("contest:an unmeasured neighbour leaves the attribution untouched",
           on["UNMEASURED"]["tier"] == "attributed_to_neighbour",
           on["UNMEASURED"]["tier"])
    record("contest:a control is never a finding, new tier included",
           on_ctrl["NTC"]["tier"] == "control_background", on_ctrl["NTC"]["tier"])

    # Off, every one of them is the single behaviour that existed before.
    record("contest:flag off leaves all four attributed",
           all(off[k]["tier"] == "attributed_to_neighbour"
               for k in ("PEDESTAL", "SIBLING", "WEAK", "UNMEASURED")),
           "; ".join(f'{k}={off[k]["tier"]}' for k in
                     ("PEDESTAL", "SIBLING", "WEAK", "UNMEASURED")))


# ---------------------------------------------------------------------------
# 17. Gates must be REACHABLE, not merely strict
# ---------------------------------------------------------------------------

def test_gate_reachability(tmp: Path) -> None:
    """A perfect finding at each plausible read length must reach the top tier.

    THIS IS THE TEST THE SUITE DID NOT HAVE, AND ITS ABSENCE COST THE PROJECT ITS WHOLE
    TOP TIER ON REAL DATA. `min_alignment_length` was a flat 100 while a 2x101 library
    trims and soft-clips to a mean of 98-99, so on real 2x100 reads nothing could ever
    be confirmed: 81 of 82 rows in the Zymo isolate batch failed it, and every
    confirmed_in_silico call in nine runs was on 150 bp SIMULATED reads. The suite could
    not see it because its one mean_alignment_length fixture is set to exactly 100.0 -
    the single value that passes.

    Every other test here asserts that a gate REJECTS bad evidence. This one asserts
    that gates ACCEPT good evidence, which is the half that was missing, and it is
    parameterised over read length precisely because that is the axis the simulated
    benchmark holds fixed at 150.
    """
    # 50 IS IN THIS LIST BECAUSE THE FIRST VERSION OF THIS TEST STOPPED AT 75, AND THE
    # FIX IT WAS WRITTEN TO GUARD REPRODUCED THE DEFECT AT 2x50. A flat illumina floor
    # of 60 is unreachable for a ~49 bp mean aligned length exactly as a flat 100 was
    # for 2x101 - the same bug, one read length down, in the test meant to prevent it.
    # That is why the floor is now a fraction of the measured read length.
    for rl, plat in [(50, "illumina"), (75, "illumina"), (101, "illumina"),
                     (151, "illumina"), (251, "illumina"), (4000, "ont")]:
        aln = int(rl * 0.98)          # trimming + soft-clipping, as observed on real data
        ident = 0.9950 if plat == "illumina" else 0.9200
        cov = tmp / f"reach_{plat}_{rl}.tsv"
        write_cov(cov, [
            cov_row(sample=f"R{rl}", role="sample", taxid="632",
                    taxon="Yersinia pestis", group="bacteria", platform=plat,
                    reference_bp=4650000, unique_reads=900000,
                    reads_per_million=90000.0, breadth=0.9500, mean_depth=29.0,
                    depth_cv=0.40, window_bp=1000, windows_expected=4650,
                    windows_occupied=4400, window_occupancy=0.9462,
                    mean_identity=ident, mean_alignment_length=aln,
                    best_neighbour_taxid="633", best_neighbour_reads=10,
                    neighbour_ratio="90000.0",
                    screen_assigned_reads=950000,
                    mapped_fraction_of_assigned="0.9500",
                    total_reads_screened=5000000, mapq_min=10),
        ])
        # The pipeline measures read length with seqkit_stats and now passes it to the
        # scorer; the fixture has to supply it too, or the test silently exercises the
        # unknown-read-length fallback instead of the path that ships.
        stats = tmp / f"stats_{plat}_{rl}.tsv"
        stats.write_text("file\tavg_len\tsample\trole\tplatform\n"
                         f"x.fastq.gz\t{rl}.0\tR{rl}\tsample\t{plat}\n")
        got = score(tmp / f"reach{rl}{plat}", cov,
                    extra=["--read-stats", stats])[f"R{rl}"]
        record(f"reachable:{plat} {rl}bp perfect evidence confirms",
               got["tier"] == "confirmed_in_silico",
               f'{got["tier"]}: {got["reasons"][:150]}')


def test_marker_block_required_to_confirm(tmp: Path) -> None:
    """A taxon with no marker rule block cannot reach confirmed_in_silico (D17a).

    31 of 90 watchlist entries have no block. For those the veto/vaccine/confirm logic
    never executes, so the finding used to be stamped confirmed on alignment geometry
    alone with an EMPTY marker_verdict and no caveat. Dormant only while nothing could
    reach the top tier; the first real-data confirmation unblocked by the reachability
    fix above was Salmonella enterica, notifiable, which had no block.

    IT HAS ONE NOW (D40), which is why this test no longer uses it: invA + ttrC + ssaV
    are real curated markers and the real zymo_rl100 library confirms on them. The
    stand-in is Orientia tsutsugamushi - notifiable, bacterial, and as of 2026-10-20 the
    only notifiable bacterium left with no block, so if it ever gains one this test
    fails loudly rather than quietly asserting the cap on a taxon that is no longer
    capped. That is exactly how this check surfaced: it was written around Salmonella
    and went red the moment Salmonella was curated.
    """
    import json
    cov = tmp / "d17a.tsv"
    def perfect(sample, taxid, taxon):
        return cov_row(sample=sample, role="sample", taxid=taxid, taxon=taxon,
                       group="bacteria", platform="illumina", reference_bp=4650000,
                       unique_reads=900000, reads_per_million=90000.0, breadth=0.9500,
                       mean_depth=29.0, depth_cv=0.40, window_bp=1000,
                       windows_expected=4650, windows_occupied=4400,
                       window_occupancy=0.9462, mean_identity=0.9950,
                       mean_alignment_length=148, best_neighbour_taxid="633",
                       best_neighbour_reads=10, neighbour_ratio="90000.0",
                       screen_assigned_reads=950000,
                       mapped_fraction_of_assigned="0.9500",
                       total_reads_screened=5000000, mapq_min=10)
    # 784 = Orientia tsutsugamushi, which really has no block in the shipped pack.
    write_cov(cov, [perfect("NOBLOCK", "784", "Orientia tsutsugamushi")])
    got = score(tmp / "d17a", cov, markers=_empty_markers(tmp))["NOBLOCK"]
    record("d17a:no marker block caps a perfect finding at supported",
           got["tier"] == "supported", got["tier"])
    record("d17a:the cap says why, and how to lift it",
           "no marker rule block" in got["reasons"]
           and "no_markers_available" in got["reasons"], got["reasons"][-160:])

    # THE REMEDY THE CAP PRINTS HAS TO BE A REMEDY THAT WORKS.
    #
    # It did not. The message above told the operator to declare
    # `no_markers_available: true`, and load_rules()'s allow-list did not contain that
    # key, so following the instruction failed the whole run with
    # "unknown rule key 'no_markers_available'". Nothing caught it because THIS test
    # built a rules file to exercise the declaration and then never passed it to
    # score() - the helper always used the shipped pack - and its assertion was
    # satisfied by the word appearing in the CAP MESSAGE itself. A test can name the
    # behaviour it means to pin and still pin nothing.
    rules_ok = tmp / "rules_declared.yaml"
    rules_ok.write_text('"Orientia tsutsugamushi":\n'
                        "    no_markers_available: true\n")
    got2 = score(tmp / "d17a_declared", cov, markers=_empty_markers(tmp),
                 rules=rules_ok)["NOBLOCK"]
    record("d17a:the declared-gap escape hatch parses at all",
           got2["tier"] != "", "an unknown-key RuleError kills the run")
    record("d17a:a declared gap lifts the cap",
           got2["tier"] == "confirmed_in_silico", got2["tier"])
    # Lifting the cap must not remove the disclosure: a declared gap and a
    # marker-confirmed call would otherwise look identical on the page.
    record("d17a:and the call still says the marker tier never ran",
           "DECLARED UNAVAILABLE" in got2["reasons"]
           and "alignment geometry alone" in got2["reasons"],
           got2["reasons"][-200:])

    # `bool("false")` IS TRUE, WHICH IS THE WHOLE REASON THIS IS PARSED STRICTLY.
    # Every other value in this file is a string, so a flag read the same way would
    # make `no_markers_available: false` grant exactly the permission it denies.
    rules_no = tmp / "rules_declined.yaml"
    rules_no.write_text('"Orientia tsutsugamushi":\n'
                        "    no_markers_available: false\n")
    got3 = score(tmp / "d17a_declined", cov, markers=_empty_markers(tmp),
                 rules=rules_no)["NOBLOCK"]
    record("d17a:`no_markers_available: false` does NOT lift the cap",
           got3["tier"] == "supported", got3["tier"])

    rules_bad = tmp / "rules_typo.yaml"
    rules_bad.write_text('"Orientia tsutsugamushi":\n'
                         "    no_markers_available: flase\n")
    pbad = run([sys.executable, BIN / "score_evidence.py", "--coverage", cov,
                "--watchlist", PACK / "watchlist.resolved.tsv",
                "--pack-yaml", PACK / "pack.yaml",
                "--marker-evidence", _empty_markers(tmp),
                "--marker-rules", rules_bad,
                "--out", tmp / "d17a_typo.tsv"])
    record("d17a:a mistyped flag is refused, not read as permission",
           pbad.returncode != 0 and "must be a boolean" in pbad.stderr,
           pbad.stderr[-200:])


def _empty_markers(tmp: Path) -> Path:
    """A marker-state file with no observations: gating ACTIVE, evidence absent.

    Needed because the cap must apply when marker gating is running and the taxon has
    no rule - not when gating is switched off altogether, where capping every
    confirmation would be a surprise rather than a safeguard.
    """
    f = tmp / "no_markers.tsv"
    f.write_text("sample\tmarker_id\tstate\ttaxon\ttier\n"
                 "NOBLOCK\tplaceholder\tUNKNOWN\tYersinia pestis\tdna\n")
    return f


def test_control_subtraction_audit(tmp: Path) -> None:
    """Control subtraction must be readable off a table, not regexed out of prose.

    D21. Across six published batches control subtraction moved 108 rows and the only
    record of any of it was the `reasons` column. Three things were invisible:

      - WHICH control decided a suppression, in machine-readable form;
      - that 11 suppressed rows carried MORE of the organism than the control did and
        were held down purely by the 5x margin, not by the contamination model;
      - that a control can carry a select agent at a level that suppresses nothing
        TODAY while standing as a floor the moment a sample lands just under it.
        lod_v3_threearm's NTC holds B. anthracis at 32.58 RPM and nothing said so.

    Each world gets its OWN coverage file and its own score() call. Mixing worlds into
    one coverage file poisons run-global `control_rpm` for every other row in it, which
    is the fixture trap that cost a day during the D12 work.
    """
    shape = dict(taxid="1639", taxon="Listeria monocytogenes", group="bacteria",
                 reference_bp=2900000, reference_contigs=1, breadth=0.9100,
                 mean_depth=12.0, depth_cv=0.40, window_bp=1000, windows_expected=2900,
                 windows_occupied=2700, window_occupancy=0.9310, mean_identity=0.9900,
                 mean_alignment_length=95, neighbour_ratio="inf",
                 mapped_fraction_of_assigned="0.9900", total_reads_screened=1000000,
                 mapq_min=10)

    def world(name, ctrl_rpm, sample_rpm):
        cov = tmp / f"audit_{name}.tsv"
        write_cov(cov, [
            cov_row(sample="NTC", role="control", unique_reads=int(ctrl_rpm * 10),
                    reads_per_million=ctrl_rpm,
                    screen_assigned_reads=int(ctrl_rpm * 10) + 5, **shape),
            cov_row(sample="SPECIMEN", role="sample", unique_reads=int(sample_rpm * 10),
                    reads_per_million=sample_rpm,
                    screen_assigned_reads=int(sample_rpm * 10) + 5, **shape),
        ])
        out = tmp / f"audit_{name}_out.tsv"
        proc = run([sys.executable, BIN / "score_evidence.py", "--coverage", cov,
                    "--watchlist", PACK / "watchlist.resolved.tsv",
                    "--pack-yaml", PACK / "pack.yaml", "--control-ratio", "5.0",
                    "--out", out])
        with out.open() as fh:
            rows = {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}
        return rows, proc.stderr

    # (a) control 100 RPM, sample 200 RPM. 200 < 5 x 100, so suppressed - while the
    #     sample holds TWICE what the control does. The margin decided this, not the
    #     contamination.
    rows, err = world("exceeds", 100.0, 200.0)
    got = rows.get("SPECIMEN", {})
    record("ctrlaudit:a sample above its control is still suppressed by the margin",
           got.get("tier") == "suppressed_control", f"got {got.get('tier')}")
    record("ctrlaudit:and the finding says so in words",
           "carries MORE of the organism" in got.get("reasons", ""),
           got.get("reasons", "")[:120])
    record("ctrlaudit:control_subtraction marks it exceeds_control",
           got.get("control_subtraction", "").startswith("exceeds_control:NTC@"),
           got.get("control_subtraction", "<missing>"))
    record("ctrlaudit:the audit table counts the exceed-control row",
           "EXCEED-CONTROL" in err, err[-300:])

    # (b) control 100, sample 50. Suppressed, and NOT an exceed case - the plain reading
    #     of "does not exceed 5x" is true here, so the warning must NOT appear.
    rows, err = world("plain", 100.0, 50.0)
    got = rows.get("SPECIMEN", {})
    record("ctrlaudit:a sample below its control is plainly suppressed",
           got.get("tier") == "suppressed_control"
           and got.get("control_subtraction", "").startswith("suppressed:NTC@"),
           f"{got.get('tier')} / {got.get('control_subtraction')}")
    record("ctrlaudit:and carries no exceed-control warning",
           "carries MORE of the organism" not in got.get("reasons", ""), "")

    # (c) control 10, sample 5000. Clears the margin. `cleared:` is NOT the same as an
    #     empty column: it records that a control held this taxon and lost.
    rows, err = world("cleared", 10.0, 5000.0)
    got = rows.get("SPECIMEN", {})
    record("ctrlaudit:a sample clearing the margin records cleared:, not empty",
           got.get("tier") != "suppressed_control"
           and got.get("control_subtraction", "").startswith("cleared:NTC@"),
           f"{got.get('tier')} / {got.get('control_subtraction')}")

    # (d) THE LATENT CASE, and the reason this test exists. The control carries the
    #     organism but suppresses nothing, because the one sample row for it is below
    #     the read floor and never reaches the suppression branch. A table printing only
    #     realized suppressions would report this control as clean while it holds a
    #     watchlist organism and imposes a standing floor.
    cov = tmp / "audit_latent.tsv"
    write_cov(cov, [
        cov_row(sample="NTC", role="control", unique_reads=300,
                reads_per_million=30.0, screen_assigned_reads=305, **shape),
        cov_row(sample="SPECIMEN", role="sample", unique_reads=20,
                reads_per_million=2.0, screen_assigned_reads=25, **shape),
    ])
    out = tmp / "audit_latent_out.tsv"
    proc = run([sys.executable, BIN / "score_evidence.py", "--coverage", cov,
                "--watchlist", PACK / "watchlist.resolved.tsv",
                "--pack-yaml", PACK / "pack.yaml", "--control-ratio", "5.0",
                "--out", out])
    record("ctrlaudit:a control that suppresses nothing is still reported, as LATENT",
           "LATENT" in proc.stderr and "Listeria monocytogenes" in proc.stderr,
           proc.stderr[-300:])
    record("ctrlaudit:the implied floor is printed, not left to be multiplied out",
           "150.00" in proc.stderr, "want 5.0 x 30.00 = 150.00 in the audit table")
    record("ctrlaudit:a row pre-empted by an earlier gate is named as such",
           "pre-empted by an earlier gate" in proc.stderr, proc.stderr[-300:])
    # A row the control never got to vote on must NOT look like a row no control
    # touched. These two were the same empty string, which made the column unable to
    # express the case the whole audit exists for.
    with (tmp / "audit_latent_out.tsv").open() as fh:
        latent = {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}
    record("ctrlaudit:a pre-empted row records not_evaluated, not empty",
           latent.get("SPECIMEN", {}).get("control_subtraction", "").startswith(
               "not_evaluated:NTC@"),
           latent.get("SPECIMEN", {}).get("control_subtraction", "<missing>"))

    # (e) THE CONTROL NAME MUST BE THE CONTROL'S NAME. The first version of the table
    #     padded to a constant and sliced [-37:], truncating from the LEFT with no
    #     marker, so d12fine's 38-40 character control names all printed as identifiers
    #     that do not exist - `od_gram_negative_fine_illumina_r1_NTC` greps to nothing.
    #     In the one table built to name the control responsible, that reads as a data
    #     error rather than a display one.
    long_name = "lod_gram_negative_fine_illumina_r1_NTC_and_then_some_more"
    cov = tmp / "audit_longname.tsv"
    write_cov(cov, [
        cov_row(sample=long_name, role="control", unique_reads=1000,
                reads_per_million=100.0, screen_assigned_reads=1005, **shape),
        cov_row(sample="SPECIMEN", role="sample", unique_reads=2000,
                reads_per_million=200.0, screen_assigned_reads=2005, **shape),
    ])
    proc = run([sys.executable, BIN / "score_evidence.py", "--coverage", cov,
                "--watchlist", PACK / "watchlist.resolved.tsv",
                "--pack-yaml", PACK / "pack.yaml", "--control-ratio", "5.0",
                "--out", tmp / "audit_longname_out.tsv"])
    # Left truncation drops the START of the name, so the discriminator is whether the
    # beginning survived - not whether some fragment appears. ("od_gram_negative" is a
    # substring of this very fixture's name, which is how the first version of this
    # check managed to fail against correct code.)
    record("ctrlaudit:a long control name keeps its beginning, not its end",
           long_name[:24] in proc.stderr, proc.stderr[-300:])
    record("ctrlaudit:and is visibly elided when it does not fit",
           long_name in proc.stderr or "\u2026" in proc.stderr, proc.stderr[-300:])

    # (f) A taxon with a control but no sample row at all is UNEXERCISED, not LATENT.
    #     Calling it LATENT buries the dangerous case - rows present, every one
    #     pre-empted - among harmless ones.
    cov = tmp / "audit_unexercised.tsv"
    write_cov(cov, [
        cov_row(sample="NTC", role="control", unique_reads=300,
                reads_per_million=30.0, screen_assigned_reads=305, **shape),
    ])
    proc = run([sys.executable, BIN / "score_evidence.py", "--coverage", cov,
                "--watchlist", PACK / "watchlist.resolved.tsv",
                "--pack-yaml", PACK / "pack.yaml", "--control-ratio", "5.0",
                "--out", tmp / "audit_unexercised_out.tsv"])
    record("ctrlaudit:a control with no sample row to act on reads UNEXERCISED",
           "UNEXERCISED" in proc.stderr and "LATENT" not in proc.stderr,
           proc.stderr[-300:])


def test_cross_scenario_control_bleed(tmp: Path) -> None:
    """D16a's mechanism, preserved as a test rather than as a batch on a disk.

    `simulate_benchmark.py` now refuses to BUILD a batch mixing a hopped control with
    other scenarios, so the batch that would demonstrate this is unconstructible by
    design. That makes this fixture the only standing record of what is at stake: a
    control belonging to one scenario suppressing genuine detections in another, because
    `score_evidence.py` keys `control_rpm` on taxid alone across the whole run.

    The numbers are the measured ones - injecting a 1,000 RPM control row into the real
    zymo_rl100 coverage set turns a 1,816.7 RPM *supported* call into
    `suppressed_control`. If anyone relaxes the hop guard, or revisits scoping
    `control_rpm` by scenario, this is the behaviour that moves.
    """
    shape = dict(taxid="1639", taxon="Listeria monocytogenes", group="bacteria",
                 reference_bp=2900000, reference_contigs=1, breadth=0.9100,
                 mean_depth=12.0, depth_cv=0.40, window_bp=1000, windows_expected=2900,
                 windows_occupied=2700, window_occupancy=0.9310, mean_identity=0.9900,
                 mean_alignment_length=95, neighbour_ratio="inf",
                 mapped_fraction_of_assigned="0.9900", total_reads_screened=1000000,
                 mapq_min=10)
    cov = tmp / "crossbleed.tsv"
    write_cov(cov, [
        # scenB's deliberately contaminated NTC. Correct in isolation.
        cov_row(sample="scenB_NTC", role="control", unique_reads=10000,
                reads_per_million=1000.0, screen_assigned_reads=10005, **shape),
        # scenA's samples, a different scenario entirely. Both genuinely positive.
        cov_row(sample="scenA_s1", role="sample", unique_reads=18167,
                reads_per_million=1816.7, screen_assigned_reads=18200, **shape),
        cov_row(sample="scenA_s2", role="sample", unique_reads=3555,
                reads_per_million=355.5, screen_assigned_reads=3600, **shape),
    ])
    got = score(tmp, cov)
    for s_id in ("scenA_s1", "scenA_s2"):
        tier = got.get(s_id, {}).get("tier", "<missing>")
        record(f"crossbleed:{s_id} is suppressed by another scenario's control",
               tier == "suppressed_control", f"got {tier}")
    record("crossbleed:the 1816 RPM row is flagged as exceeding its suppressor",
           got.get("scenA_s1", {}).get("control_subtraction", "").startswith(
               "exceeds_control:scenB_NTC@"),
           got.get("scenA_s1", {}).get("control_subtraction", "<missing>"))
    record("crossbleed:and the suppressing control is named, not just its RPM",
           "scenB_NTC" in got.get("scenA_s2", {}).get("reasons", ""),
           got.get("scenA_s2", {}).get("reasons", "")[:120])


def test_index_hop_batch_guard(tmp: Path) -> None:
    """D16a: a deliberately contaminated control cannot share a batch.

    The scenario is correct alone and incoherent in company - its NTC holds real target
    reads by design, and `score_evidence.py` subtracts controls run-globally. Refused at
    the batch level rather than scoped in the scorer, because that run-global reach is
    exactly what D14's matrix blank depends on: re-scoring lod_v4 without its controls
    turns 15 T. gondii and 15 E. coli rows back into false positives.

    Asserted in BOTH directions. A guard that refuses everything passes the first check
    and is the likelier mistake, since the hop scenario must stay runnable alone.
    """
    sim = SCRIPTS / "simulate_benchmark.py"
    head = ("defaults:\n  replicates: 1\n  platforms: [illumina]\n  host: human\n"
            "  host_fraction: 0.95\n  total_reads: 100000\nscenarios:\n")
    hop = ("  - id: hop_one\n    description: hopped control\n"
           "    kind: contamination\n    spikes:\n"
           "        - taxon: Yersinia pestis\n          source: target\n"
           "          levels_rpm: [5000]\n    truth:\n"
           "        present: [Yersinia pestis]\n    controls:\n"
           "        include_ntc: true\n        index_hop_fraction: 0.001\n")

    def plain(sid):
        return (f"  - id: {sid}\n    description: an ordinary series\n"
                f"    kind: lod_series\n    spikes:\n"
                f"        - taxon: Bacillus anthracis\n          source: target\n"
                f"          levels_rpm: [100]\n    truth:\n"
                f"        present: [Bacillus anthracis]\n")

    def attempt(text, extra, tag):
        f = tmp / f"hopguard_{tag}.yaml"
        f.write_text(text)
        proc = run([sys.executable, sim, "--scenarios", f, "--batch", "guardtest",
                    "--dry-run"] + list(extra))
        return proc.stdout + proc.stderr

    MSG = "cannot share a batch"
    mixed = attempt(head + hop + plain("plain_two"), [], "mixed")
    record("hopguard:a hop scenario sharing a batch is refused", MSG in mixed,
           mixed[-200:])
    record("hopguard:the refusal names the escape hatch",
           "--scenario hop_one" in mixed, mixed[-200:])
    record("hopguard:the hop scenario alone is NOT refused",
           MSG not in attempt(head + hop + plain("plain_two"),
                              ["--scenario", "hop_one"], "alone"), "")
    record("hopguard:two ordinary scenarios are NOT refused",
           MSG not in attempt(head + plain("plain_two") + plain("plain_three"),
                              [], "two_plain"), "")
    record("hopguard:a hop scenario alone in its own file is NOT refused",
           MSG not in attempt(head + hop, [], "solo"), "")


def test_scorer_defaults_match_shipped_config(tmp: Path) -> None:
    """An argparse default that disagrees with nextflow.config is a trap for re-scoring.

    D18 set `min_neighbour_uniformity = 0` in nextflow.config after measuring the gate on
    real reads. The argparse default stayed 0.50, so every manual re-score of a published
    batch that forgot the flag silently reproduced PRE-D18 behaviour. That happened during
    the D16a review (2026-10-13) and produced 7 tier changes where the shipped
    configuration produces 5 - a discrepancy that read as an error in docs/decisions.md
    which was not in fact there.

    Only thresholds that exist on BOTH sides are compared; a param with no CLI
    counterpart is not a mismatch. The final check fails if the pairs stop matching the
    files at all, so this cannot silently degrade into testing nothing.
    """
    import re as _re
    cfg = (REPO / "nextflow.config").read_text()
    src = (REPO / "bin" / "score_evidence.py").read_text()
    checked = 0
    for param, flag in [("min_neighbour_uniformity", "--min-neighbour-uniformity"),
                        ("pedestal_max_uniformity", "--pedestal-max-uniformity"),
                        ("control_ratio", "--control-ratio")]:
        m_cfg = _re.search(rf"^\s*{param}\s*=\s*([0-9.]+)", cfg, _re.M)
        m_arg = _re.search(rf'add_argument\(\s*"{flag}".*?default=([0-9.]+)',
                           src, _re.S)
        if not m_cfg or not m_arg:
            continue
        checked += 1
        record(f"defaults:{param} argparse default matches nextflow.config",
               float(m_cfg.group(1)) == float(m_arg.group(1)),
               f"config says {m_cfg.group(1)}, argparse says {m_arg.group(1)}")
    record("defaults:the config/argparse pairs still resolve", checked >= 2,
           f"only {checked} pair(s) matched both files")

    # THE LIST ABOVE MUST NOT SILENTLY STOP BEING THE WHOLE LIST. A new numeric param
    # that exists on both sides and is not named above would escape the check entirely,
    # which is how the min_neighbour_uniformity mismatch survived in the first place.
    # Derive the pairs from the files and require the hand-written list to match.
    cfg_params = set(_re.findall(r"^\s*([a-z_]+)\s*=\s*[0-9.]+\s*$", cfg, _re.M))
    cli_flags = set(_re.findall(r'add_argument\(\s*"--([a-z-]+)"[^)]*?default=[0-9.]',
                                src, _re.S))
    both = {f.replace("-", "_") for f in cli_flags} & cfg_params
    named = {"min_neighbour_uniformity", "pedestal_max_uniformity", "control_ratio"}
    record("defaults:every numeric param on both sides is actually checked",
           both <= named, f"unchecked pair(s): {sorted(both - named)}")


def test_protein_tier_honours_catalogue(tmp: Path) -> None:
    """The protein tier must use the floor the curator wrote, not a blanket discount.

    D23. `aa_min_id = max(0.50, min_id_catalogue - 0.15)` subtracted 15 identity points
    from every marker regardless of curation. The catalogue already encodes the
    distinction the slack was trying to express: housekeeping genes that cross-match
    across genera carry floors of 0.95-0.98, and the divergent-virus block carries
    0.78-0.90 ON PURPOSE, with a note saying the HMM tier is what carries detection
    there. The blanket subtraction overwrote exactly that curation.

    Both directions are asserted from one catalogue, because a change that merely
    tightened everything would pass the first assertion alone - and tightening the
    filovirus markers is the specific harm, since they are the evidence D15 needs.
    """
    cat = tmp / "markers_prot.tsv"
    cat.write_text(
        "# test catalogue\n"
        "marker_id\ttaxon\tgene\taliases\tmarker_class\tevidence_type\tlocus\t"
        "min_coverage\tmin_identity\tcopy_expectation\tquery_strategy\tsource_db\tnotes\n"
        # A conserved housekeeping gene, floored high by the curator precisely because
        # it cross-matches. Observed 0.877 is below 0.98 and above the old 0.83.
        "mtb_rpoB\tMycobacterium tuberculosis\trpoB\t\trequired\tprotein\tchr\t"
        "0.60\t0.98\tsingle\texact\ttest\thousekeeping, cross-matches\n"
        # A deliberately low-floored divergent-virus marker. Observed 0.880 clears 0.85.
        "ebola_NP\tZaire ebolavirus\tNP\t\trequired\tprotein\tgenome\t"
        "0.60\t0.85\tsingle\texact\ttest\tdivergent virus, low floor on purpose\n")

    # m8: query target fident alnlen mismatch gapopen qstart qend tstart tend evalue bits
    hits = tmp / "prot.m8"
    hits.write_text(
        "r1\tmarker|mtb_rpoB|X|rpoB\t0.877\t100\t0\t0\t1\t100\t1\t100\t1e-40\t200.0\n"
        "r2\tmarker|ebola_NP|Y|NP\t0.880\t100\t0\t0\t1\t100\t1\t100\t1e-40\t200.0\n")
    faa = tmp / "prot.faa"
    faa.write_text(">marker|mtb_rpoB|X|rpoB\n" + "M" * 100 + "\n"
                   ">marker|ebola_NP|Y|NP\n" + "M" * 100 + "\n")

    # A MANIFEST IS REQUIRED FOR THIS TEST TO MEAN ANYTHING. Without one, a marker
    # that fails every tier is reported UNKNOWN rather than ABSENT - deliberately, so
    # "never built" is not mistaken for "absent from sample". A fixture with no
    # manifest therefore cannot tell a tightened floor from an unbuilt marker, and
    # would pass whatever the floor did.
    man = tmp / "markers_prot.manifest.tsv"
    man.write_text(
        "marker_id\ttaxon\tgene\tmarker_class\tevidence_type\tquery_strategy\tstatus\n"
        "mtb_rpoB\tMycobacterium tuberculosis\trpoB\trequired\tprotein\texact\tbuilt\n"
        "ebola_NP\tZaire ebolavirus\tNP\trequired\tprotein\texact\tbuilt\n")

    def states(extra):
        out = tmp / f"ms_{abs(hash(str(extra)))}.tsv"
        run([sys.executable, BIN / "marker_evidence.py", "--markers", cat,
             "--manifest", man,
             "--aa-hits", hits, "--aa-fasta", faa, "--sample", "S1",
             "--out", out] + list(extra))
        with out.open() as fh:
            return {r["marker_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    got = states([])
    record("protfloor:a housekeeping gene below its catalogue floor is ABSENT",
           got.get("mtb_rpoB", {}).get("state") == "ABSENT",
           f"got {got.get('mtb_rpoB', {}).get('state')} "
           f"(aa_id 0.877 vs catalogue 0.98)")
    record("protfloor:a divergent-virus marker above its own low floor is PRESENT",
           got.get("ebola_NP", {}).get("state") == "PRESENT",
           f"got {got.get('ebola_NP', {}).get('state')} "
           f"(aa_id 0.880 vs catalogue 0.85)")

    # The old behaviour is still reachable, declared rather than implicit.
    old = states(["--protein-identity-slack", "0.15"])
    record("protfloor:the old blanket slack is still reachable when asked for",
           old.get("mtb_rpoB", {}).get("state") == "PRESENT",
           f"got {old.get('mtb_rpoB', {}).get('state')} at slack 0.15")

    # ONT keeps a small protein allowance, and it must NOT compound with the much
    # larger nucleotide relaxation - the bug the file already records for min_id.
    ont = states(["--platform", "ont", "--ont-identity-slack", "0.20"])
    record("protfloor:ONT protein slack does not compound with the nucleotide one",
           ont.get("mtb_rpoB", {}).get("state") == "ABSENT",
           f"got {ont.get('mtb_rpoB', {}).get('state')}; 0.877 must stay below "
           f"0.98-0.05, not 0.98-0.20")

    # aa_length and hmm_qlen exist so the tiers can be re-examined from the artifact.
    record("protfloor:aa_length is emitted, so aa_coverage has a denominator",
           (got.get("ebola_NP", {}).get("aa_length") or "") not in ("", "0"),
           f"aa_length={got.get('ebola_NP', {}).get('aa_length')!r}")


def test_marker_expectation_cli(tmp: Path) -> None:
    """The second expectation axis, end to end through the command line.

    The evaluator's own self-test covers the grading rules against hand-built
    dicts. This covers the seam those cannot reach: --marker-states on the
    command line, load_marker_states reading a real marker_evidence.py file, and
    the exit code that gates CI.

    The axis exists because a findings row is only ever written for a watchlist
    taxon. divergent_clostridium_bont spikes C. butyricum, which is deliberately
    a near_neighbour and not on the watchlist, so its only expectation could be
    `max_tier: marker_veto` on C. botulinum - an upper bound, which "no finding"
    satisfies. The scenario was green on every run while the toxin it exists to
    test went unreported.
    """
    sc = tmp / "mk_scenarios.yaml"
    sc.write_text(
        "schema_version: 1\n"
        "scenarios:\n"
        "    - id: tox\n"
        "      description: a toxin in an untargeted neighbour\n"
        "      kind: specificity\n"
        "      spikes:\n"
        "          - taxon: Clostridium butyricum\n"
        "            source: neighbour\n"
        "            levels_rpm: [300]\n"
        "      truth:\n"
        "          present: [Clostridium butyricum]\n"
        "          absent: [Clostridium botulinum]\n"
        "      expect:\n"
        "          - taxon: Clostridium botulinum\n"
        "            max_tier: marker_veto\n"
        "          - marker: cbut_bont\n"
        "            state: PRESENT\n"
        "            taxon: Clostridium butyricum\n")

    man = tmp / "mk_manifest.tsv"
    cols = ("sample\tscenario_id\tkind\treplicate\tplatform\tlevel_rpm\trole\t"
            "truth_present\ttruth_absent\ttruth_vaccine\n")
    man.write_text(cols + "".join(
        f"t{i}\ttox\tspecificity\t{i}\tillumina\t300\tsample\t"
        f"Clostridium butyricum\tClostridium botulinum\t\n" for i in (1, 2, 3)))

    # No findings row for C. butyricum is CORRECT - it is not on the watchlist.
    # That is precisely why the old expectation could not fail.
    fnd = tmp / "mk_findings.tsv"
    fnd.write_text("sample\trole\tplatform\ttier\ttaxid\ttaxon\n")

    def states(path, value):
        path.write_text(
            "sample\tmarker_id\ttaxon\tgene\tmarker_class\tstate\ttier_used\t"
            "tiers_not_built\treason\n"
            + "".join(f"t{i}\tcbut_bont\tClostridium butyricum\tbont\t"
                      f"discriminating\t{value}\tdna,protein\thmm\tx\n"
                      for i in (1, 2, 3)))
        return path

    def run_eval(extra, tag):
        out = tmp / f"mk_eval_{tag}"
        proc = run([sys.executable, BIN / "evaluate_benchmark.py",
                    "--scenarios", sc, "--manifest", man, "--findings", fnd,
                    "--out-dir", out] + list(extra))
        rows = []
        f = out / "expectations.tsv"
        if f.is_file():
            with f.open() as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
        return proc, {r["check"]: r for r in rows}

    # The toxin is absent from the sample: the scenario must now go red.
    absent = states(tmp / "mk_absent.tsv", "ABSENT")
    proc, got = run_eval(["--marker-states", absent], "absent")
    record("markerexp:a marker expectation can FAIL a scenario findings cannot",
           got.get("marker_state[cbut_bont]", {}).get("result") == "FAIL",
           f"got {got.get('marker_state[cbut_bont]', {}).get('result')!r}")
    record("markerexp:the old taxon expectation still passes on no finding",
           got.get("max_tier", {}).get("result") == "PASS",
           f"got {got.get('max_tier', {}).get('result')!r} - if this flipped, the "
           f"new axis changed the old one")
    record("markerexp:a failing marker expectation gates CI",
           proc.returncode != 0, f"exit {proc.returncode}")
    record("markerexp:every replicate is named in the detail",
           all(f"t{i}" in (got.get("marker_state[cbut_bont]", {}).get("detail") or "")
               for i in (1, 2, 3)),
           f"detail={got.get('marker_state[cbut_bont]', {}).get('detail')!r}")

    # The same scenario against a sample that does carry the toxin.
    present = states(tmp / "mk_present.tsv", "PRESENT")
    proc_ok, got_ok = run_eval(["--marker-states", present], "present")
    record("markerexp:the toxin reported on every replicate PASSes",
           got_ok.get("marker_state[cbut_bont]", {}).get("result") == "PASS",
           f"got {got_ok.get('marker_state[cbut_bont]', {}).get('result')!r}")
    record("markerexp:and the run then exits clean",
           proc_ok.returncode == 0, f"exit {proc_ok.returncode}")

    # An unbuilt marker is a database gap. It must not read as a detection
    # failure, for the same reason marker_evidence.py refuses to call it ABSENT.
    unk = states(tmp / "mk_unknown.tsv", "UNKNOWN")
    proc_u, got_u = run_eval(["--marker-states", unk], "unknown")
    record("markerexp:an UNKNOWN marker SKIPs rather than failing the scenario",
           got_u.get("marker_state[cbut_bont]", {}).get("result") == "SKIP",
           f"got {got_u.get('marker_state[cbut_bont]', {}).get('result')!r}")
    record("markerexp:the unbuilt tier is named, so the SKIP is explicable",
           "hmm" in (got_u.get("marker_state[cbut_bont]", {}).get("detail") or ""),
           f"detail={got_u.get('marker_state[cbut_bont]', {}).get('detail')!r}")

    # Omitting the flag entirely must not turn every marker expectation red:
    # that would make --marker-states mandatory for every caller.
    proc_n, got_n = run_eval([], "none")
    record("markerexp:without --marker-states the expectation SKIPs, not FAILs",
           got_n.get("marker_state[cbut_bont]", {}).get("result") == "SKIP",
           f"got {got_n.get('marker_state[cbut_bont]', {}).get('result')!r}")
    record("markerexp:and a run with no marker states still exits clean",
           proc_n.returncode == 0, f"exit {proc_n.returncode}")


def test_marker_tier_gap(tmp: Path) -> None:
    """A tier that was never built must not be reported as having been searched.

    `tiers_supplied` answers "did the run get an input file for this tier". It was
    being used as if it answered "does this marker have anything in that file",
    and `tried` was appended unconditionally inside each tier's branch. A marker
    whose HMM was never built therefore emitted `tier_used = dna,hmm,protein` and
    a reason reading "searched via dna,hmm,protein; no tier met thresholds".

    On the shipped pack that was 13 markers, and the missing tier is almost always
    the HMM one - cbut_bont, ebola_NP, marburg_GP, nipah_L, rvfv_L - every one a
    divergence-catching marker on a high-consequence agent, which is exactly the
    set markers.tsv argues the HMM tier exists for. orthopox_ATI and mpox_C3L are
    missing the NUCLEOTIDE tier instead, which is worse: a Variola marker whose
    reason claimed a nucleotide search that never happened.

    The verdict must NOT move. This is a reporting defect; an unbuilt tier has no
    hits either way, so a fix that changed a state would be changing something
    else as well.
    """
    cat = tmp / "gap_markers.tsv"
    cat.write_text(
        "marker_id\ttaxon\tgene\taliases\tmarker_class\tevidence_type\tlocus\t"
        "min_coverage\tmin_identity\tcopy_expectation\tquery_strategy\tsource_db\tnotes\n"
        # Declares an HMM tier it does not have, and is not in the sample.
        "gap_absent\tTestus absentus\tg1\t\trequired\tprotein+hmm\tchr\t"
        "0.60\t0.80\tsingle\texact\ttest\tdeclares hmm, has none\n"
        # Same declaration, same gap, but the protein tier finds it.
        "gap_present\tTestus praesens\tg2\t\trequired\tprotein+hmm\tchr\t"
        "0.60\t0.80\tsingle\texact\ttest\tdeclares hmm, has none, is present\n"
        # Fully built: the negative control that keeps the fix from being a
        # blanket "never report hmm".
        "full_absent\tTestus completus\tg3\t\trequired\tprotein+hmm\tchr\t"
        "0.60\t0.80\tsingle\texact\ttest\teverything built\n"
        # Every declared tier is missing, so nothing was searched at all.
        "none_built\tTestus nihil\tg4\t\trequired\thmm\tchr\t"
        "0.60\t0.80\tsingle\texact\ttest\tonly tier is unbuilt\n")

    hits = tmp / "gap.m8"
    hits.write_text(
        "r1\tmarker|gap_present|Y|g2\t0.95\t100\t0\t0\t1\t100\t1\t100\t1e-40\t200.0\n")
    faa = tmp / "gap.faa"
    faa.write_text(">marker|gap_absent|X|g1\n" + "M" * 100 + "\n"
                   ">marker|gap_present|Y|g2\n" + "M" * 100 + "\n"
                   ">marker|full_absent|Z|g3\n" + "M" * 100 + "\n")
    # Supplied but empty: the HMM tier's INPUT exists for this run, which is what
    # made the old code believe every marker had been searched by it.
    domtbl = tmp / "gap.domtbl"
    domtbl.write_text("# hmmsearch domtbl, no hits\n")

    man = tmp / "gap.manifest.tsv"
    man.write_text(
        "marker_id\ttaxon\tgene\tmarker_class\tevidence_type\tquery_strategy\t"
        "status\tn_nt_seqs\tn_aa_seqs\thas_hmm\n"
        "gap_absent\tTestus absentus\tg1\trequired\tprotein+hmm\texact\tpartial\t0\t2\tno\n"
        "gap_present\tTestus praesens\tg2\trequired\tprotein+hmm\texact\tpartial\t0\t2\tno\n"
        "full_absent\tTestus completus\tg3\trequired\tprotein+hmm\texact\tbuilt\t0\t2\tyes\n"
        "none_built\tTestus nihil\tg4\trequired\thmm\texact\tpartial\t0\t0\tno\n")

    # The same manifest WITHOUT the per-tier columns. A pack built before they
    # existed must keep its old behaviour rather than reading as all-gaps.
    old_man = tmp / "gap.oldmanifest.tsv"
    old_man.write_text(
        "marker_id\ttaxon\tgene\tmarker_class\tevidence_type\tquery_strategy\tstatus\n"
        "gap_absent\tTestus absentus\tg1\trequired\tprotein+hmm\texact\tpartial\n"
        "gap_present\tTestus praesens\tg2\trequired\tprotein+hmm\texact\tpartial\n"
        "full_absent\tTestus completus\tg3\trequired\tprotein+hmm\texact\tbuilt\n"
        "none_built\tTestus nihil\tg4\trequired\thmm\texact\tpartial\n")

    def states(manifest, tag):
        out = tmp / f"gapstates_{tag}.tsv"
        run([sys.executable, BIN / "marker_evidence.py", "--markers", cat,
             "--manifest", manifest, "--aa-hits", hits, "--aa-fasta", faa,
             "--hmm-hits", domtbl, "--sample", "S1", "--out", out])
        with out.open() as fh:
            return {r["marker_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    got = states(man, "new")

    record("tiergap:an unbuilt tier is not listed in tier_used",
           got.get("gap_absent", {}).get("tier_used") == "protein",
           f"got {got.get('gap_absent', {}).get('tier_used')!r}, want 'protein'")
    record("tiergap:the unbuilt tier is named in its own column",
           got.get("gap_absent", {}).get("tiers_not_built") == "hmm",
           f"got {got.get('gap_absent', {}).get('tiers_not_built')!r}")
    record("tiergap:the ABSENT reason says the tier was not searched",
           "NOT SEARCHED" in (got.get("gap_absent", {}).get("reason") or "")
           and "hmm" in (got.get("gap_absent", {}).get("reason") or ""),
           f"reason={got.get('gap_absent', {}).get('reason')!r}")
    record("tiergap:the verdict does not move - reporting fix only",
           got.get("gap_absent", {}).get("state") == "ABSENT",
           f"got {got.get('gap_absent', {}).get('state')}")

    record("tiergap:a PRESENT row still carries the caveat",
           got.get("gap_present", {}).get("state") == "PRESENT"
           and "NOT SEARCHED" in (got.get("gap_present", {}).get("reason") or ""),
           f"state={got.get('gap_present', {}).get('state')}, "
           f"reason={got.get('gap_present', {}).get('reason')!r}")

    record("tiergap:a fully built marker still reports every tier",
           got.get("full_absent", {}).get("tier_used") == "hmm,protein"
           and not got.get("full_absent", {}).get("tiers_not_built"),
           f"tier_used={got.get('full_absent', {}).get('tier_used')!r}, "
           f"not_built={got.get('full_absent', {}).get('tiers_not_built')!r}")

    # Nothing was searched, so absence cannot be asserted. ABSENT here would be a
    # database gap masquerading as evidence against the organism - the same class
    # of error as the mtb_RD1 vaccine misfire this module already records.
    record("tiergap:a marker with no built tier is UNKNOWN, not ABSENT",
           got.get("none_built", {}).get("state") == "UNKNOWN",
           f"got {got.get('none_built', {}).get('state')}")
    record("tiergap:and says the tier is missing rather than the input",
           "not built" in (got.get("none_built", {}).get("reason") or ""),
           f"reason={got.get('none_built', {}).get('reason')!r}")

    old = states(old_man, "old")
    record("tiergap:a manifest without the columns keeps the old behaviour",
           old.get("gap_absent", {}).get("tier_used") == "hmm,protein"
           and not old.get("gap_absent", {}).get("tiers_not_built"),
           f"tier_used={old.get('gap_absent', {}).get('tier_used')!r}, "
           f"not_built={old.get('gap_absent', {}).get('tiers_not_built')!r}")
    # The one place the fix DOES move a verdict, pinned in both directions so the
    # change is deliberate rather than incidental. Without the columns the old code
    # searched a tier that does not exist, found nothing, and called that ABSENT -
    # a database gap reported as evidence against the organism. With them it is
    # UNKNOWN, which cannot promote or veto.
    record("tiergap:the old manifest asserts an absence it never searched for",
           old.get("none_built", {}).get("state") == "ABSENT",
           f"got {old.get('none_built', {}).get('state')}")
    record("tiergap:the per-tier columns turn that false absence into UNKNOWN",
           got.get("none_built", {}).get("state") == "UNKNOWN",
           f"got {got.get('none_built', {}).get('state')}")


def test_protein_coverage_per_sequence(tmp: Path) -> None:
    """Protein coverage is per allele; the summed denominator made it unreachable.

    `fasta_lengths()` returned marker_id -> SUMMED length and `aa_evidence()`
    divided the union of target intervals by it. A marker holding 8 alternative
    representations of one locus therefore capped at 0.125 however well any
    single allele was covered, against catalogue floors of 0.50-0.60.

    Measured on the shipped pack when this was found: 107 of 125 protein markers
    could not reach their own `min_coverage` on a PERFECT full-length hit -
    ebola_NP, nipah_L, marburg_NP, lassa_NP and most of the other
    divergent-agent markers among them. The failure is not a missing call: the
    row writes a confident ABSENT, which feeds `absent()` as TRUE and can fire a
    veto_if or vaccine_if rule.

    This is the same defect nt_evidence() documents fixing ("BREADTH IS PER
    SEQUENCE, NOT OVER THEIR SUMMED LENGTH ... structurally unreachable at ANY
    depth"), which was never carried to the protein tier. Third place it has
    been found; the suite had nothing pinning it, and passed throughout.
    """
    cat = tmp / "cov_markers.tsv"
    cat.write_text(
        "marker_id\ttaxon\tgene\taliases\tmarker_class\tevidence_type\tlocus\t"
        "min_coverage\tmin_identity\tcopy_expectation\tquery_strategy\tsource_db\tnotes\n"
        "multi\tTestus multiplex\tg1\t\trequired\tprotein\tchr\t"
        "0.60\t0.80\tsingle\texact\ttest\tfour alleles of one locus\n"
        "single\tTestus solus\tg2\t\trequired\tprotein\tchr\t"
        "0.60\t0.80\tsingle\texact\ttest\tone allele, the control\n")
    man = tmp / "cov.manifest.tsv"
    man.write_text(
        "marker_id\ttaxon\tgene\tmarker_class\tevidence_type\tquery_strategy\t"
        "status\tn_nt_seqs\tn_aa_seqs\thas_hmm\n"
        "multi\tTestus multiplex\tg1\trequired\tprotein\texact\tbuilt\t0\t4\tno\n"
        "single\tTestus solus\tg2\trequired\tprotein\texact\tbuilt\t0\t1\tno\n")

    # Four alleles of 100 aa each: summed 400, best 100.
    faa = tmp / "cov.faa"
    faa.write_text(
        "".join(f">marker|multi|A{i} allele {i} [Testus multiplex]\n" + "M" * 100 + "\n" for i in range(4))
        + ">marker|single|B0 only [Testus solus]\n" + "M" * 100 + "\n")

    def states(hits_text, tag):
        hits = tmp / f"cov_{tag}.m8"
        hits.write_text(hits_text)
        out = tmp / f"covstates_{tag}.tsv"
        run([sys.executable, BIN / "marker_evidence.py", "--markers", cat,
             "--manifest", man, "--aa-hits", hits, "--aa-fasta", faa,
             "--sample", "S1", "--out", out])
        with out.open() as fh:
            return {r["marker_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    # One allele covered END TO END. Under the summed denominator this is 100/400
    # = 0.25 and ABSENT; per sequence it is 1.0 and PRESENT.
    full = states("r\tmarker|multi|A2\t0.95\t100\t0\t0\t1\t100\t1\t100\t1e-9\t200.0\n"
                  "r\tmarker|single|B0\t0.95\t100\t0\t0\t1\t100\t1\t100\t1e-9\t200.0\n",
                  "full")
    record("aacov:a full-length hit on ONE allele is PRESENT",
           full.get("multi", {}).get("state") == "PRESENT",
           f"got {full.get('multi', {}).get('state')!r}, "
           f"coverage {full.get('multi', {}).get('aa_coverage')!r}")
    record("aacov:coverage is 1.0, not 1/n_alleles",
           full.get("multi", {}).get("aa_coverage") == "1.0000",
           f"got {full.get('multi', {}).get('aa_coverage')!r}, want '1.0000'")
    record("aacov:the denominator is one allele, not their sum",
           full.get("multi", {}).get("aa_length") == "100",
           f"got {full.get('multi', {}).get('aa_length')!r}, want '100' (sum is 400)")
    record("aacov:the chosen allele is named",
           full.get("multi", {}).get("aa_best_ref") == "marker|multi|A2",
           f"got {full.get('multi', {}).get('aa_best_ref')!r}")
    record("aacov:a single-allele marker is unaffected",
           full.get("single", {}).get("state") == "PRESENT",
           f"got {full.get('single', {}).get('state')!r}")

    # THE GATE MUST STILL REJECT. A fix that makes everything reachable by
    # shrinking the denominator to nothing would pass the checks above.
    part = states("r\tmarker|multi|A2\t0.95\t20\t0\t0\t1\t20\t1\t20\t1e-9\t200.0\n",
                  "partial")
    record("aacov:a 20% hit on that same allele is still ABSENT",
           part.get("multi", {}).get("state") == "ABSENT",
           f"got {part.get('multi', {}).get('state')!r}, "
           f"coverage {part.get('multi', {}).get('aa_coverage')!r}")

    # Coverage may not be accumulated ACROSS alleles to clear the floor: 30% of
    # each of four alleles is not 120% of the locus.
    spread = states("".join(
        f"r\tmarker|multi|A{i}\t0.95\t30\t0\t0\t1\t30\t1\t30\t1e-9\t200.0\n"
        for i in range(4)), "spread")
    record("aacov:partial hits on four alleles do not sum to coverage",
           spread.get("multi", {}).get("state") == "ABSENT",
           f"got {spread.get('multi', {}).get('state')!r}, "
           f"coverage {spread.get('multi', {}).get('aa_coverage')!r}")

    # A target absent from the FASTA has no length, so its coverage silently
    # computes as 0. nt_evidence() warns on exactly this; the protein tier did not.
    hits = tmp / "cov_unknown.m8"
    hits.write_text("r\tmarker|multi|ZZZ\t0.95\t100\t0\t0\t1\t100\t1\t100\t1e-9\t200.0\n")
    out = tmp / "covstates_unknown.tsv"
    pr = run([sys.executable, BIN / "marker_evidence.py", "--markers", cat,
              "--manifest", man, "--aa-hits", hits, "--aa-fasta", faa,
              "--sample", "S1", "--out", out])
    record("aacov:a target missing from the FASTA is warned about, not silent",
           "absent from the marker FASTA" in pr.stderr,
           f"stderr tail: {pr.stderr[-200:]!r}")


def test_marker_catalogue_watchlist_join(tmp: Path) -> None:
    """Three hand-maintained files agree on literal strings; something must check.

    markers.tsv.taxon, marker_rules.yaml's keys and watchlist.resolved.tsv.name name the
    same organisms. A drift in any one degrades SILENTLY - a rule block whose key stops
    matching simply stops gating, which is indistinguishable from a taxon that has no
    rules, which is D17a's permissive path.
    """
    proc = run([sys.executable, BIN / "marker_rules.py",
                "--rules", PACK / "marker_rules.yaml",
                "--markers", PACK / "markers.tsv",
                "--watchlist", PACK / "watchlist.resolved.tsv", "--validate"])
    record("markerjoin:the shipped pack validates", proc.returncode == 0,
           proc.stderr[-300:])
    record("markerjoin:a rule block with no marker of its own is reported",
           "Brucella suis" in proc.stderr, proc.stderr[-300:])
    record("markerjoin:a watchlist taxon with markers but no rule block is reported",
           "Monkeypox virus" in proc.stderr, proc.stderr[-300:])

    # A rule block keyed on a name no watchlist row carries can never fire. That must
    # be an ERROR, not a note - it is silent permissiveness.
    bad = tmp / "bad_rules.yaml"
    bad.write_text('"Not A Real Organism":\n'
                   '    confirm_requires: mtb_rpoB\n'
                   '    rationale: test\n')
    proc = run([sys.executable, BIN / "marker_rules.py", "--rules", bad,
                "--markers", PACK / "markers.tsv",
                "--watchlist", PACK / "watchlist.resolved.tsv", "--validate"])
    record("markerjoin:a rule block that is not a watchlist name FAILS validation",
           proc.returncode != 0 and "can never fire" in proc.stderr,
           f"rc={proc.returncode} {proc.stderr[-200:]}")


def test_lead_join_seam(tmp: Path) -> None:
    """resolve_watchlist writes the descendants file; evaluate_benchmark reads it.

    Two files, one hand-written column contract, and a rename on either side degrades
    SILENTLY: load_covers() would return an empty mapping, every group lead would read
    as an unlisted false positive, and every species it covered would go back to being
    a flat miss - with no error anywhere. The evaluator's own self-test checks the join
    logic against a hand-built dict; this checks that a real file feeds it.
    """
    sys.path.insert(0, str(BIN))
    import evaluate_benchmark as eb

    # The contract, asserted on a file written by hand so it fails loudly on a
    # rename rather than skipping when no taxdump is around.
    d = tmp / "desc_min.tsv"
    d.write_text("watchlist_taxid\twatchlist_name\tdescendant_taxid\tdescendant_name\n"
                 "11266\tFiloviridae\t565995\tBundibugyo ebolavirus\n"
                 "11266\tFiloviridae\t186538\tZaire ebolavirus\n"
                 # A retired id with no current scientific name. It stays in the file
                 # for parse_screen, which matches taxids, and must be dropped here
                 # rather than inserting an empty name that matches nothing.
                 "11266\tFiloviridae\t999999\t\n")
    covers = eb.load_covers(d)
    record("leadjoin:load_covers reads the file resolve_watchlist writes",
           covers.get("Filoviridae") == {"Bundibugyo ebolavirus", "Zaire ebolavirus"},
           f"got {dict(covers)}")
    record("leadjoin:a descendant with no scientific name is dropped",
           "" not in covers.get("Filoviridae", set()), f"got {dict(covers)}")
    record("leadjoin:a missing descendants file is not an error",
           eb.load_covers(None) == {} and eb.load_covers(tmp / "nope.tsv") == {},
           "load_covers must return empty, not raise - --descendants is optional")

    # Round trip against the real taxdump: the writer's header must still parse.
    real = tmp / "desc_real.tsv"
    proc = run([sys.executable, BIN / "resolve_watchlist.py",
                "--watchlist", PACK / "watchlist.tsv",
                "--taxdump", os.environ.get(
                    "BTF_TAXDUMP", "/media/8tb/bioThreatFinder/taxdump/2026-09-29"),
                "--out", tmp / "wl_rt.tsv", "--descendants-out", real])
    if proc.returncode != 0 or not real.exists():
        skip("leadjoin:the real writer round-trips", "no usable taxdump (set BTF_TAXDUMP)")
        skip("leadjoin:descendant names are real scientific names", "no usable taxdump")
        return
    rt = eb.load_covers(real)
    # THE NAME THE TRUTH USES IS OFTEN NOT THE SCIENTIFIC NAME. NCBI files
    # "Bundibugyo ebolavirus" as an equivalent name of "Bundibugyo virus", and the
    # benchmark manifest, like the published outbreak lists it was built from, uses
    # the equivalent. Emitting scientific names alone made the join miss every one of
    # them silently: the group read as an unlisted false positive and the species it
    # covered stayed a flat miss. This is the string that caught it.
    filo = rt.get("Filoviridae", set())
    if filo:
        record("leadjoin:an equivalent name joins, not just the scientific name",
               "Bundibugyo ebolavirus" in filo and "Bundibugyo virus" in filo,
               f"both must be present; got {sorted(n for n in filo if 'undibugyo' in n)}")
    else:
        skip("leadjoin:an equivalent name joins, not just the scientific name",
             "no Filoviridae entry in the shipped watchlist")
    record("leadjoin:the real writer round-trips",
           bool(rt) and all(isinstance(v, set) and v for v in rt.values()),
           f"{len(rt)} entries covered")
    # Casing matters: the join is against watchlist and manifest strings, which carry
    # NCBI casing. A lowercased name would match nothing and fail silently.
    mm = rt.get("Marburg marburgvirus", set())
    record("leadjoin:descendant names are real scientific names",
           any(n[:1].isupper() for n in mm) if mm else False,
           f"sample of {len(mm)}: {sorted(mm)[:3]}")


WATCH_COLS = ["taxid", "name", "resolved_rank", "group", "category",
              "min_unique_reads", "min_breadth", "near_neighbour_taxids",
              "n_descendants", "notes", "source", "resolution_note"]


def _watch(path: Path, rows: list[dict]) -> Path:
    out = [ "\t".join(WATCH_COLS) ]
    for r in rows:
        unknown = sorted(set(r) - set(WATCH_COLS))
        if unknown:
            raise KeyError(f"_watch() got column(s) not in WATCH_COLS: {unknown}")
        out.append("\t".join(str(r.get(c, "")) for c in WATCH_COLS))
    path.write_text("\n".join(out) + "\n")
    return path


def test_supra_specific_entries_are_leads(tmp: Path) -> None:
    """A family cannot be confirmed as an organism, however good the geometry.

    A watchlist entry resolved above species exists so a divergent member with no
    reference of its own still surfaces - Kraken2 rolls Bundibugyo's reads up through
    Filoviridae and every species-keyed row stays silent. What that supports is "a
    member of this group is here", which is a lead. Scoring it as a species call would
    put a family name in the `confirmed_in_silico` column.

    The cap is explicit rather than left to breadth arithmetic, because a family
    bundle's breadth divides the covered span by the summed ref_bp of every genome
    sharing the taxid - so whether it clears `min_occupancy` depends on how many
    genomes the family happens to have, and adding one more would change the verdict.
    """
    wl = _watch(tmp / "ranks.tsv", [
        # Identical category and thresholds; `resolved_rank` is the only variable.
        dict(taxid="11266", name="Filoviridae", resolved_rank="family", group="virus",
             category="high_consequence", min_unique_reads=20, min_breadth=0.30,
             n_descendants=103, source="test"),
        dict(taxid="632", name="Yersinia pestis", resolved_rank="species",
             group="bacteria", category="high_consequence", min_unique_reads=100,
             min_breadth=0.01, near_neighbour_taxids="633;630", n_descendants=287,
             source="test"),
        # 18 of the 89 shipped entries resolve to NCBI rank "no rank" and every one
        # names a real organism. A `!= "species"` test would cap Variola.
        dict(taxid="10255", name="Variola virus", resolved_rank="no rank",
             group="virus", category="high_consequence", min_unique_reads=20,
             min_breadth=0.30, near_neighbour_taxids="10245;10243;10244;28873;12643",
             n_descendants=4, source="test"),
    ])

    def viral(sample, taxid, taxon, **kw):
        d = dict(sample=sample, role="sample", taxid=taxid, taxon=taxon, group="virus",
                 reference_bp=19000, unique_reads=400, reads_per_million=80.0,
                 breadth=0.9200, mean_depth=4.80, depth_cv=0.90, window_bp=190,
                 windows_expected=100, windows_occupied=92, window_occupancy=0.9200,
                 mean_identity=0.9960, neighbour_ratio="inf",
                 screen_assigned_reads=420, mapped_fraction_of_assigned="0.9523",
                 total_reads_screened=5000000, mapq_min=10)
        d.update(kw)
        return cov_row(**d)

    cov = tmp / "leadcap.tsv"
    write_cov(cov, [
        viral("FAMILY", "11266", "Filoviridae"),
        viral("NORANK", "10255", "Variola virus", reference_bp=185000,
              window_bp=1850, windows_expected=100),
        # Same bacterial geometry that reaches confirmed_in_silico in test_tiering.
        cov_row(sample="SPECIES", role="sample", taxid="632", taxon="Yersinia pestis",
                group="bacteria", reference_bp=4650000, unique_reads=820,
                reads_per_million=164.0, breadth=0.3400, mean_depth=0.88, depth_cv=1.15,
                window_bp=1000, windows_expected=4650, windows_occupied=2790,
                window_occupancy=0.6000, mean_identity=0.9950,
                best_neighbour_taxid="633", best_neighbour_reads=12,
                neighbour_ratio="68.33", screen_assigned_reads=1000,
                mapped_fraction_of_assigned="0.8200", total_reads_screened=5000000,
                mapq_min=10),
    ])
    f = score(tmp, cov, extra=["--watchlist", str(wl)])

    fam = f.get("FAMILY", {})
    record("leadcap:a group-rank entry is capped at group_lead",
           fam.get("tier") == "group_lead",
           f"got {fam.get('tier')}, geometry was breadth 0.92 occupancy 0.92")
    record("leadcap:the cap names the rank and the tier it replaced",
           "CAPPED AT GROUP_LEAD" in fam.get("reasons", "")
           and "family" in fam.get("reasons", ""),
           fam.get("reasons", "")[-200:])

    sp = f.get("SPECIES", {})
    record("leadcap:identical evidence at species rank is NOT capped",
           sp.get("tier") == "confirmed_in_silico", f"got {sp.get('tier')}")
    nr = f.get("NORANK", {})
    record("leadcap:NCBI rank 'no rank' is NOT capped",
           nr.get("tier") in ("supported", "confirmed_in_silico"),
           f"got {nr.get('tier')} - Variola and 17 others resolve to 'no rank'")

    # BUNDLE DILUTION. Measured on the real run, not invented: the Filoviridae entry
    # fired on all three Bundibugyo replicates at 4,017 RPM and 5.4x depth with its
    # member covered end to end, and `breadth` read 0.0858 - which is 18,940 / 220,597,
    # the member's share of a 12-genome bundle. Under the 0.20 floor that row landed in
    # `screened_only`, the bucket for classifier noise, so the entry bought nothing an
    # analyst would see. best_reference_breadth is the measure that answers the question
    # a group entry actually asks.
    cov3 = tmp / "leadcap_bundle.tsv"
    write_cov(cov3, [
        viral("DILUTED", "11266", "Filoviridae", reference_bp=220597,
              unique_reads=7984, reads_per_million=4016.64,
              breadth=0.085813, window_occupancy=0.084071, windows_occupied=8,
              mean_depth=5.4289, depth_cv=3.3154,
              best_reference_breadth=0.998, best_reference_occupancy=0.990,
              best_reference="NC_014373.1"),
    ])
    h = score(tmp, cov3, extra=["--watchlist", str(wl)])
    dil = h.get("DILUTED", {})
    record("leadcap:the substitution is recorded in the reasons",
           "GROUP GEOMETRY" in dil.get("reasons", "")
           and "NC_014373.1" in dil.get("reasons", ""),
           "the published row keeps the bundle breadth, so a reader comparing it to "
           f"the floor needs to be told: {dil.get('reasons','')[:160]}")
    record("leadcap:a group row is judged on its best reference, not the bundle",
           dil.get("tier") == "group_lead",
           f"got {dil.get('tier')}; bundle breadth 0.0858 is below the 0.20 floor, "
           f"best-reference breadth 0.998 is not")

    # The substitution is for group entries ONLY. A species whose bundle happens to
    # carry several assemblies must keep being judged on the whole thing, or a target
    # covered in one of its three genomes would read as fully covered.
    cov4 = tmp / "leadcap_species_bundle.tsv"
    write_cov(cov4, [
        cov_row(sample="SPBUNDLE", role="sample", taxid="632", taxon="Yersinia pestis",
                group="bacteria", reference_bp=13950000, unique_reads=820,
                reads_per_million=164.0, breadth=0.0050, mean_depth=0.30, depth_cv=4.10,
                window_bp=1000, windows_expected=13950, windows_occupied=70,
                window_occupancy=0.0050, mean_identity=0.9950,
                best_neighbour_taxid="633", best_neighbour_reads=12,
                neighbour_ratio="68.33", screen_assigned_reads=1000,
                mapped_fraction_of_assigned="0.8200", total_reads_screened=5000000,
                mapq_min=10, best_reference_breadth=0.980,
                best_reference_occupancy=0.975,
                best_reference="NZ_CP009973.1"),
    ])
    i = score(tmp, cov4, extra=["--watchlist", str(wl)])
    record("leadcap:a SPECIES row is still judged on the whole bundle",
           i.get("SPBUNDLE", {}).get("tier") not in ("supported", "confirmed_in_silico"),
           f"got {i.get('SPBUNDLE', {}).get('tier')}; 0.005 bundle breadth must not be "
           f"rescued by one well-covered assembly")

    # THE CAP MUST ONLY LOWER. `marker_veto` and `suppressed_control` outrank
    # `screened_only` in TIER_ORDER but assert the OPPOSITE of presence; lowering them
    # would turn a recorded negative into an open lead, and a control row reclassified
    # to `control_background` would be promoted into a finding.
    cov2 = tmp / "leadcap_floor.tsv"
    write_cov(cov2, [
        viral("SUPPRESSED", "11266", "Filoviridae", reads_per_million=30.0,
              unique_reads=150),
        viral("NTC", "11266", "Filoviridae", role="control", reads_per_million=22.0,
              unique_reads=110, breadth=0.4000, windows_occupied=40,
              window_occupancy=0.4000, total_reads_screened=2500000),
    ])
    g = score(tmp, cov2, extra=["--watchlist", str(wl)])
    record("leadcap:a suppressed group row is not lifted to a lead",
           g.get("SUPPRESSED", {}).get("tier") == "suppressed_control",
           f"got {g.get('SUPPRESSED', {}).get('tier')}")
    record("leadcap:a family row in a control stays control_background",
           g.get("NTC", {}).get("tier") == "control_background",
           f"got {g.get('NTC', {}).get('tier')}")

    # The shipped pack holds nothing supra-specific, so the cap changes no finding
    # today. If that stops being true this assertion is the notice.
    with (PACK / "watchlist.resolved.tsv").open() as fh:
        ranks = {r["name"]: r["resolved_rank"]
                 for r in csv.DictReader(fh, delimiter="\t")}
    sys.path.insert(0, str(BIN))
    import score_evidence
    supra = sorted(n for n, r in ranks.items()
                   if r.strip().lower() in score_evidence.SUPRA_SPECIFIC_RANKS)
    # Filoviridae was added deliberately (D25) and is the only one. Pinning the exact
    # list, not just "non-empty": a group entry caps every finding under it at a lead,
    # so one arriving by accident - a watchlist name that happens to resolve to a
    # genus - would silently stop a species from ever being called, and this is the
    # only place that would say so.
    record("leadcap:exactly the intended supra-specific entries are present",
           supra == ["Filoviridae"], f"capped entries: {supra}")


def _crashed_sylph_flag(tmp: Path, k2: Path, wl: Path, desc: Path) -> str:
    """The flag sylph gets when its output file is absent entirely (it crashed),
    as opposed to present and empty (it found nothing)."""
    out = tmp / "summary_sylph_dead.tsv"
    run([sys.executable, BIN / "parse_screen.py", "--kraken2-report", k2,
         "--expected-arms", "kraken2,sylph", "--watchlist", wl,
         "--descendants", desc, "--sample", "S4",
         "--out", tmp / "cand4.tsv", "--out-summary", out])
    if not out.exists():
        return "<no summary>"
    with out.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["classifier"] == "sylph":
                return r.get("flag", "")
    return "<no sylph row>"


def test_report_screen_coverage(tmp: Path) -> None:
    """The report must say what stage 1 FAILED to search, not only what it found.

    Nothing in this suite exercised make_report.py before, and the thing it was
    silently not doing is the reason: screen_summary.tsv was published to the
    results directory and read by nothing. The two facts that live only in it are
    the two a reader most needs in order not to over-read a clean page -

      * an arm that was ENABLED and produced nothing. METABULI, KRAKENUNIQ and
        SYLPH run under `errorStrategy 'ignore'` so that one dead arm cannot sink
        the run; the cost is that its failure is indistinguishable from a clean
        negative. For a taxon something else also raised, the missing arm shows up
        in findings.tsv's `classifiers` column - but a taxon ONLY the dead arm
        would have raised leaves no trace anywhere, and a report that does not say
        an arm died is claiming a search that did not happen.

      * the unclassified fraction. Perseus (Bioinformatics 2026) measured that
        sequence from a family absent from the database is not left unclassified
        but CONFIDENTLY ASSIGNED to a known taxon, approaching 50% on long reads.
        The novel agent this tool exists to catch is then not missed; it is
        mislabelled as background. The unclassified fraction is the only signal
        that survives that, so it has to be on the page.
    """
    findings = tmp / "findings.tsv"
    findings.write_text(
        "sample\trole\tplatform\ttier\ttaxid\ttaxon\tunique_reads\treasons\n"
        "S1\tsample\tillumina\tconfirmed_in_silico\t632\tYersinia pestis\t5000\tgates met\n")

    cols = ("sample\trole\tplatform\tclassifier\tstatus\ttotal_reads\t"
            "unclassified_reads\tunclassified_fraction\tflag\tnote\n")
    screen = tmp / "screen_summary.tsv"
    screen.write_text(cols + "".join([
        "S1\tsample\tillumina\tkraken2\tok\t1000000\t550000\t0.550000\tHIGH_UNCLASSIFIED\t\n",
        "S1\tsample\tillumina\tmetabuli\tno_output\t\t\t\tARM_PRODUCED_NOTHING\tenabled, no report\n",
        "S1\tsample\tillumina\tbracken\tno_output\t\t\t\t\troutine on a clean negative\n",
        "S1\tsample\tillumina\tsylph\tok\t\t\t\t\t\n",
        "S1\tsample\tillumina\tkrakenuniq\tnot_run\t\t\t\t\tdisabled for this run\n",
    ]))

    html, summ = tmp / "r.html", tmp / "s.tsv"
    p = run([sys.executable, BIN / "make_report.py", "--findings", findings,
             "--screen-summary", screen, "--pack", "human",
             "--html", html, "--summary", summ])
    if not html.exists():
        record("report:renders with a screen summary", False, p.stderr[-300:])
        return
    doc = html.read_text()
    record("report:renders with a screen summary", True)

    record("report:a dead arm is bannered at the top, not buried in a table",
           "produced nothing: metabuli" in doc
           and "NOT a negative result" in doc
           and doc.index("produced nothing: metabuli") < doc.index("Stage-1 screen"),
           "the banner is what an analyst reads before deciding the page is clean")
    record("report:the unclassified fraction reaches the page",
           "0.550000" in doc and "HIGH_UNCLASSIFIED" in doc)
    record("report:a high unclassified fraction is bannered",
           "unclassified (up to 55%)" in doc, doc[doc.find("high fraction"):][:120])
    record("report:every arm appears in the coverage table",
           all(a in doc for a in ("kraken2", "metabuli", "bracken", "sylph",
                                  "krakenuniq")))
    record("report:bracken's routine absence carries no alarm pill",
           "pill crit'>ARM" in doc and doc.count("pill crit'>ARM") == 1,
           "only metabuli's failure should be pilled, not bracken's by-design exit")

    # THE RUN-LEVEL VIEW IS WHERE THE SUPPRESSION ABOVE IS PAID FOR.
    #
    # parse_screen.py deliberately does not flag a pack-restricted arm's empty
    # profile, because on one sample it is the correct negative. The cost of that is
    # a real blind spot: an empty syldb, or one truncated at build time, LOADS
    # WITHOUT ERROR and returns zero rows for every sample, which per-sample is
    # indistinguishable from a cohort that genuinely carries nothing. Across the run
    # it is not indistinguishable at all, and the report is the only place that holds
    # every sample at once. An arm that matched nothing anywhere is bannered; an arm
    # that matched nothing on some samples and something on others is the normal
    # mixed cohort and must stay silent.
    multi = tmp / "screen_multi.tsv"
    multi.write_text(cols + "".join([
        f"S{i}\tsample\tillumina\tkraken2\tok\t1000\t10\t0.010000\t\t\n"
        f"S{i}\tsample\tillumina\tsylph\tempty_report\t\t\t\t\tnegative answer\n"
        for i in (1, 2, 3)]))
    html3 = tmp / "r3.html"
    run([sys.executable, BIN / "make_report.py", "--findings", findings,
         "--screen-summary", multi, "--pack", "human",
         "--html", html3, "--summary", tmp / "s3.tsv"])
    doc3 = html3.read_text() if html3.exists() else ""
    record("report:an arm empty on EVERY sample is bannered as a database problem",
           "matched nothing anywhere in this run" in doc3
           and "sylph on all 3 sample(s)" in doc3,
           "an empty or truncated database is silent per-sample and loud per-run")
    record("report:an arm empty on every sample is not called a crash",
           "produced nothing: sylph" not in doc3,
           "it ran; the question is whether its database had anything in it")

    mixed = tmp / "screen_mixed.tsv"
    mixed.write_text(cols + "".join([
        "S1\tsample\tillumina\tsylph\tempty_report\t\t\t\t\tnegative answer\n",
        "S2\tsample\tillumina\tsylph\tok\t\t\t\t\t\n",
        "S3\tsample\tillumina\tsylph\tempty_report\t\t\t\t\tnegative answer\n",
    ]))
    html4 = tmp / "r4.html"
    run([sys.executable, BIN / "make_report.py", "--findings", findings,
         "--screen-summary", mixed, "--pack", "human",
         "--html", html4, "--summary", tmp / "s4.tsv"])
    doc4 = html4.read_text() if html4.exists() else ""
    record("report:a mixed cohort raises no database alarm",
           bool(doc4) and "matched nothing anywhere" not in doc4,
           "one `ok` row proves the database loaded; the empties are then just "
           "samples without a watchlist organism in them")

    # WITHOUT THE SUMMARY THE REPORT MUST LOOK THE SAME AS IT ALWAYS DID.
    # The file is `optional:` nowhere now, but a report rendered from an older run
    # directory has no such file, and the section must simply be absent rather than
    # render an empty table that reads as "every arm was fine".
    html2 = tmp / "r2.html"
    run([sys.executable, BIN / "make_report.py", "--findings", findings,
         "--pack", "human", "--html", html2, "--summary", tmp / "s2.tsv"])
    record("report:no summary means no section, not an empty one implying all-clear",
           html2.exists() and "Stage-1 screen coverage" not in html2.read_text())

    # A PRESENTATIONAL STEP MUST NOT BE ABLE TO END THE RUN.
    #
    # MAKE_REPORT sets no errorStrategy, so an exception here kills the pipeline
    # AFTER stage 1, stage 2, the marker tier and scoring have all finished - the
    # actual result traded for the view of it. bf9744d fixed exactly that in
    # PAVIAN_EXPORT, and wiring a new input into the LAST step of the run re-opens
    # it. A non-numeric unclassified_fraction is the cheap way in.
    bad = tmp / "screen_bad.tsv"
    bad.write_text(cols + "S1\tsample\tillumina\tkraken2\tok\tn/a\tn/a\t"
                          "not-a-number\tHIGH_UNCLASSIFIED\t\n")
    html3 = tmp / "r3.html"
    p = run([sys.executable, BIN / "make_report.py", "--findings", findings,
             "--screen-summary", bad, "--pack", "human",
             "--html", html3, "--summary", tmp / "s3.tsv"])
    record("report:a malformed screen summary degrades, it does not kill the run",
           p.returncode == 0 and html3.exists(),
           f"rc={p.returncode}: {p.stderr.strip()[-200:]}")

    # And a file that cannot be parsed at all must SAY so on the page. Silence
    # there is the same error one level up: a reader takes a missing coverage
    # table for a screen with nothing to report.
    unreadable = tmp / "screen_unreadable.tsv"
    unreadable.write_bytes(b"\xff\xfe\x00binary garbage, not a TSV\x00\xff")
    html4 = tmp / "r4.html"
    p = run([sys.executable, BIN / "make_report.py", "--findings", findings,
             "--screen-summary", unreadable, "--pack", "human",
             "--html", html4, "--summary", tmp / "s4.tsv"])
    ok = p.returncode == 0 and html4.exists()
    doc4 = html4.read_text() if html4.exists() else ""
    record("report:an unreadable screen summary still produces a report", ok,
           f"rc={p.returncode}: {p.stderr.strip()[-200:]}")
    record("report:and the page says the coverage is unknown rather than clean",
           "could not be read" in doc4 and "clean screen" in doc4,
           "a missing table reads as 'nothing to report'; the page has to say "
           "the arms are UNKNOWN, not omit the section as it does when no summary "
           "was supplied at all")



def test_samplesheet_role_authority(tmp: Path) -> None:
    """The samplesheet declares which samples are controls, so it must DRIVE subtraction.

    D39. `control_rpm` - the dict that actually suppresses anything - was built only from
    the `role` field on the pipeline's own output rows, while `n_controls`, which decides
    the QC sentence an analyst reads first, also honoured --samplesheet. Inside one run
    those agree, because both descend from the same Nextflow `meta`. They diverge on the
    RE-SCORE path this project uses routinely: published coverage TSVs plus a corrected
    sheet. Declaring a control that way silenced "NO NEGATIVE CONTROL IN RUN" and left
    control_rpm empty, selecting the "controls all clean" branch - so every finding
    asserted no taxon had been raised in any control. On zymo_rl100 that sentence went
    onto all 82 findings while the same file reported the declared control as a
    `supported` L. monocytogenes detection at 355.5 RPM.

    Each world gets its own coverage file: run-global `control_rpm` is poisoned by mixing
    worlds, which is the fixture trap recorded on test_control_subtraction_audit.
    """
    shape = dict(taxid="1639", taxon="Listeria monocytogenes", group="bacteria",
                 reference_bp=2900000, reference_contigs=1, breadth=0.9100,
                 mean_depth=12.0, depth_cv=0.40, window_bp=1000, windows_expected=2900,
                 windows_occupied=2700, window_occupancy=0.9310, mean_identity=0.9900,
                 mean_alignment_length=95, neighbour_ratio="inf",
                 mapped_fraction_of_assigned="0.9900", total_reads_screened=1000000,
                 mapq_min=10)

    def sheet(path: Path, rows: list[tuple[str, str]]) -> Path:
        path.write_text("sample,platform,fastq_1,fastq_2,role\n" + "".join(
            f"{s},illumina,/x/{s}_1.fq.gz,/x/{s}_2.fq.gz,{r}\n" for s, r in rows))
        return path

    def stats(path: Path, samples: list[str]) -> Path:
        path.write_text("file\tavg_len\tsample\trole\tplatform\n" + "".join(
            f"/x/{s}.fq.gz\t101.0\t{s}\tsample\tillumina\n" for s in samples))
        return path

    def world(name, cov_rows, *, ss=None, rs=None):
        cov = tmp / f"ssrole_{name}.tsv"
        write_cov(cov, cov_rows)
        out = tmp / f"ssrole_{name}_out.tsv"
        cmd = [sys.executable, BIN / "score_evidence.py", "--coverage", cov,
               "--watchlist", PACK / "watchlist.resolved.tsv",
               "--pack-yaml", PACK / "pack.yaml", "--control-ratio", "5.0",
               "--out", out]
        if ss:
            cmd += ["--samplesheet", ss]
        if rs:
            cmd += ["--read-stats", rs]
        proc = run(cmd)
        with out.open() as fh:
            rows = {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}
        return rows, proc.stderr

    def rows_for(ntc_role):
        return [
            cov_row(sample="NTC", role=ntc_role, unique_reads=1000,
                    reads_per_million=100.0, screen_assigned_reads=1005, **shape),
            cov_row(sample="SPECIMEN", role="sample", unique_reads=2000,
                    reads_per_million=200.0, screen_assigned_reads=2005, **shape),
        ]

    # (a) THE BUG. The rows say NTC is an ordinary sample - the shape of every re-score -
    #     while the sheet declares it a control. 200 RPM < 5 x 100, so the specimen row
    #     must be suppressed, and the NTC's own 100 RPM must become QC rather than a
    #     finding about a specimen.
    ss = sheet(tmp / "ssrole_a.csv", [("NTC", "control"), ("SPECIMEN", "sample")])
    rs = stats(tmp / "ssrole_a_stats.tsv", ["NTC", "SPECIMEN"])
    rows, err = world("declared", rows_for("sample"), ss=ss, rs=rs)
    spec = rows.get("SPECIMEN", {})
    record("ssrole:a sheet-declared control subtracts even when the rows say otherwise",
           spec.get("tier") == "suppressed_control", f"got {spec.get('tier')}")
    record("ssrole:the declared control's own row is QC, not a finding",
           rows.get("NTC", {}).get("tier") == "control_background",
           f"got {rows.get('NTC', {}).get('tier')}")
    record("ssrole:and NOTHING claims the control was clean",
           "all clean" not in spec.get("reasons", ""), spec.get("reasons", "")[:160])
    record("ssrole:the relabelling is reported, not silent",
           "samplesheet roles applied" in err and "is a CONTROL per the samplesheet" in err,
           err[-300:])

    # (b) The sheet is authority in BOTH directions. A row labelled control that the sheet
    #     calls an ordinary sample must not impose a floor on the run - otherwise a stale
    #     role on a re-scored row silently suppresses real findings.
    ss = sheet(tmp / "ssrole_b.csv", [("NTC", "sample"), ("SPECIMEN", "sample")])
    rs = stats(tmp / "ssrole_b_stats.tsv", ["NTC", "SPECIMEN"])
    rows, err = world("demoted", rows_for("control"), ss=ss, rs=rs)
    record("ssrole:a sheet-demoted control stops subtracting",
           rows.get("SPECIMEN", {}).get("tier") != "suppressed_control",
           f"got {rows.get('SPECIMEN', {}).get('tier')}")
    record("ssrole:and the run is then correctly reported as having no control",
           "NO NEGATIVE CONTROL IN RUN" in rows.get("SPECIMEN", {}).get("reasons", ""),
           rows.get("SPECIMEN", {}).get("reasons", "")[:160])

    # (c) No sheet at all: the pre-D39 behaviour must be untouched, since that is how
    #     every in-pipeline run scores.
    rows, err = world("nosheet", rows_for("control"))
    record("ssrole:with no sheet the row roles still decide, unchanged",
           rows.get("SPECIMEN", {}).get("tier") == "suppressed_control",
           f"got {rows.get('SPECIMEN', {}).get('tier')}")

    # (d) THE GOOD CASE THIS MUST NOT BREAK. A spotless NTC produces no coverage rows at
    #     all. It ran - it has a read ledger entry - so it is a clean control and the
    #     finding must say so instead of "contamination cannot be excluded".
    ss = sheet(tmp / "ssrole_d.csv", [("NTC", "control"), ("SPECIMEN", "sample")])
    rs = stats(tmp / "ssrole_d_stats.tsv", ["NTC", "SPECIMEN"])
    rows, err = world("spotless", [rows_for("sample")[1]], ss=ss, rs=rs)
    spec = rows.get("SPECIMEN", {})
    record("ssrole:a spotless declared control is still credited as clean",
           "negative control(s) in run, all clean" in spec.get("reasons", ""),
           spec.get("reasons", "")[:160])
    record("ssrole:and a clean control does not read as a missing one",
           "NO NEGATIVE CONTROL IN RUN" not in spec.get("reasons", ""), "")

    # (e) THE PHANTOM. Same empty-handed signature as (d), but the control is absent from
    #     the read ledger too - a sheet left over from another batch, or a name that never
    #     matched what the pipeline emitted. It never ran, so it excludes nothing, and
    #     crediting it would hand the analyst the strongest claim on the weakest evidence.
    ss = sheet(tmp / "ssrole_e.csv", [("NTC", "control"), ("SPECIMEN", "sample")])
    rs = stats(tmp / "ssrole_e_stats.tsv", ["SPECIMEN"])
    rows, err = world("phantom", [rows_for("sample")[1]], ss=ss, rs=rs)
    spec = rows.get("SPECIMEN", {})
    record("ssrole:a declared control that never ran is NOT credited",
           "NO NEGATIVE CONTROL IN RUN" in spec.get("reasons", ""),
           spec.get("reasons", "")[:160])
    record("ssrole:and it is named on stderr rather than dropped quietly",
           "declares 'NTC' a negative control" in err and "never ran" in err,
           err[-300:])

    # (f) Without a read ledger there is nothing to distinguish (d) from (e), so the
    #     question must be left unasked - the spotless case keeps its credit.
    ss = sheet(tmp / "ssrole_f.csv", [("NTC", "control"), ("SPECIMEN", "sample")])
    rows, err = world("noledger", [rows_for("sample")[1]], ss=ss)
    record("ssrole:with no read ledger a declared control is given the benefit of doubt",
           "negative control(s) in run, all clean"
           in rows.get("SPECIMEN", {}).get("reasons", ""),
           rows.get("SPECIMEN", {}).get("reasons", "")[:160])



# Markers whose built nucleotide reference is a DIFFERENT GENE from the one the
# catalogue declares. Each is a real defect in the shipped marker database, each needs
# a re-retrieval, and each is listed here with its cause so that the invariant below can
# guard against NEW ones without going red on the known two. D40.
#
# This list is self-retiring, not an allowlist that rots: the test asserts that every
# entry STILL mismatches. Rebuild the marker database with the fixed matcher and these
# rows come back correct, the assertion fails, and the list must be emptied. A pending
# item cannot outlive the problem it describes.
# (a) The builder's OWN matcher would not accept these gene tags for these markers -
#     the reference is a different gene. EMPTY as of the 2026-10-21 rebuild: `nipah_L`
#     was the only entry and no longer holds the wrong coding sequence.
GENE_TAG_WRONG_GENE: dict = {}
# (b) The matcher accepts the tag, but the DECLARED gene is not among the retrieved
#     ones. Not all of these are defects, so each says which it is.
GENE_TAG_NO_EXACT = {
    "mers_ORF1a": "NOT a defect: ORF1a is the 5' portion of the ORF1ab polyprotein and "
                  "NCBI annotates the CDS as ORF1ab. Correct sequence, different name. "
                  "Permanent entry.",
}
# Entries that must be FIXED, as opposed to accepted. Empty: both original entries were
# verified resolved against the rebuilt database on 2026-10-21.
GENE_TAG_MUST_FIX: set = set()


def test_marker_gene_identity(tmp: Path) -> None:
    """A marker's sequences must be the GENE THE CATALOGUE NAMES.

    D40. Nothing checked this, and three of the shipped markers failed it. The
    nucleotide and protein tiers are retrieved by separate ladders, so they can land on
    different genes while the manifest reports `status=built, nt=3, aa=3` and
    marker_evidence.py treats the marker as fully testable. `min_coverage` is then a
    fraction of the wrong reference.

    Found while adding the Salmonella block: `salm_ttrC` listed its operon partners
    ttrA/ttrB as aliases, and CDS extraction returned 3x3063 bp of ttrA as the
    nucleotide reference while the protein tier correctly fetched 340 aa of ttrC.
    """
    sys.path.insert(0, str(SCRIPTS))
    try:
        from build_marker_db import gene_match_rank as rank, gene_name_matches as matches
    except Exception as exc:
        record("markergene:build_marker_db importable", False, str(exc))
        return

    # ---- the collision that put the Nipah P gene in the L marker ----
    record("markergene:a single-letter CDS name cannot claim a word alias",
           not matches("P", ["L", "polymerase"]),
           "Nipah P matched the L marker because 'polymerase'.startswith('p')")
    record("markergene:nor claim a two-character one",
           not matches("E", ["E2", "TC-83"]), "")
    record("markergene:the gene it really is still matches",
           rank("L", ["L", "polymerase"]) == 0, f"rank={rank('L', ['L', 'polymerase'])}")

    # ---- the permissiveness that is there ON PURPOSE must survive ----
    record("markergene:an operon annotation still matches its member",
           rank("capBCADE", ["capB"]) == 1, f"rank={rank('capBCADE', ['capB'])}")
    record("markergene:case-insensitive prefix still matches",
           matches("STX2A", ["stx2"]), "")
    record("markergene:ORF1ab still answers for ORF1a",
           matches("ORF1ab", ["ORF1a", "upE"]), "")
    record("markergene:a genuine truncation of three characters still matches",
           rank("pag", ["pagA"]) == 2 and rank("inv", ["invA"]) == 2,
           f"pag={rank('pag', ['pagA'])}, inv={rank('inv', ['invA'])}")
    record("markergene:an unrelated gene is still refused",
           not matches("recA", ["capB"]), "")

    # ---- exact beats prefix, which is what caf1 needed ----
    record("markergene:an exact gene match outranks an operon-partner match",
           rank("caf1", ["caf1", "F1 capsular antigen"])
           < rank("caf1A", ["caf1", "F1 capsular antigen"]),
           "caf1 must win over caf1A, or the operon decides which genes a marker holds")

    # ---- the invariant, against the marker database that actually shipped ----
    db = Path(os.environ.get("BTF_MARKER_DB",
                             "/media/8tb/bioThreatFinder/packs/human/markers"))
    man = db / "markers.manifest.tsv"
    if not man.is_file():
        skip("markergene:built references name the declared gene",
             f"no marker database at {db} (set BTF_MARKER_DB)")
        return

    cat: dict = {}
    for line in (PACK / "markers.tsv").read_text().splitlines():
        if line.startswith("#") or line.startswith("marker_id") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) > 3:
            cat[f[0]] = (f[2], f[3])

    wrong_gene: dict = {}
    no_exact: dict = {}
    checked = 0
    with man.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            mid = r["marker_id"]
            fna = db / "nt" / f"{mid}.fna"
            if not fna.is_file():
                continue
            genes = sorted(set(re.findall(r"\[gene=([^\]]+)\]", fna.read_text())))
            if not genes:
                continue          # Route 2 took the whole record; no CDS gene to check.
            checked += 1
            declared, aliases = cat.get(mid, (r.get("gene", ""), ""))
            wanted = [declared] + [a for a in aliases.split(";") if a.strip()]
            # (a) JUDGED BY THE BUILDER'S OWN MATCHER, not by set membership. Asking a
            #     different question here than the retriever asked would make the
            #     invariant disagree with the code it is meant to constrain: 'ORF1ab'
            #     is a legitimate rank-1 answer for a wanted 'ORF1a' and must pass,
            #     while the Nipah 'P' CDS is no answer at all for 'L' and must not.
            if not any(rank(g, wanted) is not None for g in genes):
                wrong_gene[mid] = genes
            # (b) Weaker, and separate: NOTHING retrieved is an exact match to the
            #     gene or to any of its declared aliases - every sequence got in on a
            #     prefix. An alias exists precisely so a differently-annotated gene is
            #     found, so `wbeT` for a declared `wbeO1` and `MP` for `M` are right
            #     answers and must pass; what this catches is a marker that holds only
            #     its operon's OTHER genes, which is ypest_pmt1_caf1.
            elif not any(rank(g, wanted) == 0 for g in genes):
                no_exact[mid] = genes

    unexpected = {k: v for k, v in wrong_gene.items() if k not in GENE_TAG_WRONG_GENE}
    record("markergene:no built reference is a gene the matcher would reject",
           not unexpected,
           f"{checked} checked; UNEXPECTED: "
           + "; ".join(f"{k} holds {v}" for k, v in unexpected.items()))

    unexpected2 = {k: v for k, v in no_exact.items() if k not in GENE_TAG_NO_EXACT}
    record("markergene:every built reference contains the gene actually declared",
           not unexpected2,
           f"{checked} checked; UNEXPECTED: "
           + "; ".join(f"{k} holds {v}" for k, v in unexpected2.items()))

    # A pending defect must not outlive the problem it records: once the marker
    # database is rebuilt these come back correct, this fails, and the entry goes.
    seen_bad = set(wrong_gene) | set(no_exact)
    stale = sorted(GENE_TAG_MUST_FIX - seen_bad)
    record("markergene:every gene-tag defect still pending is still real",
           not stale,
           f"fixed now - remove from the pending lists: {stale}")


def test_salmonella_block(tmp: Path) -> None:
    """D40. The Salmonella block is species confirmation, and needs all three loci.

    invA and ttrC are GENUS-level - S. bongori carries both - so the species call rests
    on ssaV, the SPI-2 T3SS that S. bongori lacks. Dropping any one of the three must
    hold the finding at `supported` rather than confirming it, and an UNTESTED marker
    must hold it just as firmly as an absent one.
    """
    def salmonella(sample: str) -> str:
        # Alignment evidence taken from the real zymo_rl100 isolate library, which
        # passes every alignment gate on its own (D37). Only the markers vary.
        return cov_row(sample=sample, role="sample", platform="illumina", taxid="28901",
                       taxon="Salmonella enterica", group="bacteria",
                       reference_bp=5009666, unique_reads=24057383,
                       reads_per_million=925901.0, breadth=0.9137,
                       mean_depth=480.42, depth_cv=0.4181, window_bp=1000,
                       windows_expected=5009, windows_occupied=4697,
                       window_occupancy=0.9378, mean_identity=0.9870,
                       mean_alignment_length=99,
                       best_neighbour_taxid="562", best_neighbour_reads=1000,
                       neighbour_ratio="100.0", screen_assigned_reads=24057383,
                       mapped_fraction_of_assigned="1.0000",
                       total_reads_screened=25983000, mapq_min=10)

    cov = tmp / "salm.tsv"
    write_cov(cov, [salmonella(s) for s in ("ALL3", "NO_SSAV", "SSAV_UNTESTED")])
    mk = tmp / "salm_markers.tsv"
    mk.write_text("sample\tmarker_id\tstate\n" + "".join([
        "ALL3\tsalm_invA\tPRESENT\n", "ALL3\tsalm_ttrC\tPRESENT\n",
        "ALL3\tsalm_ssaV\tPRESENT\n",
        # Confidently not found: S. bongori would look like this.
        "NO_SSAV\tsalm_invA\tPRESENT\n", "NO_SSAV\tsalm_ttrC\tPRESENT\n",
        "NO_SSAV\tsalm_ssaV\tABSENT\n",
        # ssaV OMITTED FROM THE FILE, not ABSENT in value - never searched.
        "SSAV_UNTESTED\tsalm_invA\tPRESENT\n", "SSAV_UNTESTED\tsalm_ttrC\tPRESENT\n",
    ]))
    # THE READ LEDGER IS PART OF THE FIXTURE, not decoration. Without it the
    # alignment-length floor falls back to a flat 100 bp and a 2x101 library's 99 bp
    # mean fails it - which is exactly the gate D37 found holding this organism back,
    # and the real run cleared it only because seqkit_stats told the scorer the reads
    # were 101 bp (floor = 60% of 101 = 61). A fixture without it would assert the
    # marker block against the wrong gate and pass for the wrong reason.
    stats = tmp / "salm_stats.tsv"
    stats.write_text("file\tavg_len\tsample\trole\tplatform\n" + "".join(
        f"/x/{s}.fq.gz\t101.0\t{s}\tsample\tillumina\n"
        for s in ("ALL3", "NO_SSAV", "SSAV_UNTESTED")))
    got = score(tmp, cov, markers=mk, extra=["--read-stats", str(stats)])

    all3 = got.get(("ALL3", "Salmonella enterica"), {})
    record("salm:all three loci present confirms",
           all3.get("tier") == "confirmed_in_silico", f"got {all3.get('tier')}")
    record("salm:and the finding names the rule it satisfied",
           "salm_invA and salm_ttrC and salm_ssaV" in all3.get("reasons", ""),
           all3.get("reasons", "")[:160])

    nossav = got.get(("NO_SSAV", "Salmonella enterica"), {})
    record("salm:genus markers without the species marker do NOT confirm",
           nossav.get("tier") != "confirmed_in_silico",
           f"got {nossav.get('tier')} - invA and ttrC are both present in S. bongori")

    untested = got.get(("SSAV_UNTESTED", "Salmonella enterica"), {})
    record("salm:an UNTESTED species marker holds it just as firmly",
           untested.get("tier") != "confirmed_in_silico",
           f"got {untested.get('tier')}")



def test_marker_state_columns(tmp: Path) -> None:
    """The states file's columns must not depend on which marker sorts first.

    D40. `marker_evidence.py` took its header from `results[0]`. A record only carries
    `nt_depth`, `nt_expected_breadth`, `nt_n_seqs` and `nt_best_ref` if it reached the
    nucleotide measurement block, so a catalogue whose FIRST row is UNKNOWN - never
    built, no input for its tier - produced a header without them, and
    `extrasaction="ignore"` then dropped those four from every other row WITHOUT error,
    including rows that had measured them.

    `nt_expected_breadth` is the input to the depth gate, the rule that marker absence
    is only called when depth could have shown it. So the column that exists to stop a
    false ABSENT could vanish silently - the exact failure this file is written to
    prevent. Observed 2026-10-20 while verifying the Salmonella block against a
    three-marker database.
    """
    if not shutil.which("samtools"):
        skip("statecols:late columns survive an UNKNOWN first row", "samtools not on PATH")
        return

    ref = "marker|zz_built|ACC1|[gene=zz]"
    sam = tmp / "mc.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n"
        f"@SQ\tSN:{ref}\tLN:60\n"
        + "".join(f"r{i}\t0\t{ref}\t1\t60\t60M\t*\t0\t0\t{'A' * 60}\t"
                  f"{'I' * 60}\tNM:i:0\n" for i in range(6)))
    bam = tmp / "mc.bam"
    if run(["samtools", "sort", "-o", bam, sam]).returncode != 0:
        skip("statecols:late columns survive an UNKNOWN first row", "samtools sort failed")
        return
    run(["samtools", "index", bam])

    # aa_unbuilt sorts FIRST and is absent from the manifest, so it returns at gate 1
    # and never acquires the four nt_* keys. zz_built comes second and does.
    cat = tmp / "mc_catalogue.tsv"
    cat.write_text(
        "marker_id\ttaxon\tgene\taliases\tmarker_class\tevidence_type\tlocus\t"
        "min_coverage\tmin_identity\tcopy_expectation\tquery_strategy\tsource_db\tnotes\n"
        "aa_unbuilt\tTestus unbuiltus\taa\t-\trequired\tdna\tchromosome\t"
        "0.50\t0.90\tsingle\tgene_in_organism\ttest\tnever built on purpose\n"
        "zz_built\tTestus builtus\tzz\t-\trequired\tdna\tchromosome\t"
        "0.50\t0.90\tsingle\tgene_in_organism\ttest\tbuilt, and measured\n")
    man = tmp / "mc_manifest.tsv"
    man.write_text(
        "marker_id\ttaxon\tgene\tmarker_class\tevidence_type\tquery_strategy\t"
        "status\tn_nt_seqs\tn_aa_seqs\thas_hmm\n"
        "zz_built\tTestus builtus\tzz\trequired\tdna\tgene_in_organism\tbuilt\t1\t0\tno\n")

    out = tmp / "mc_states.tsv"
    proc = run([sys.executable, BIN / "marker_evidence.py", "--markers", cat,
                "--manifest", man, "--sample", "S1", "--nt-bam", bam,
                "--platform", "illumina", "--out", out])
    if not out.exists():
        record("statecols:late columns survive an UNKNOWN first row", False,
               proc.stderr[-300:])
        return

    with out.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        cols = list(rdr.fieldnames or [])
        rows = {r["marker_id"]: r for r in rdr}

    late = ["nt_depth", "nt_expected_breadth", "nt_n_seqs", "nt_best_ref"]
    record("statecols:late columns survive an UNKNOWN first row",
           all(c in cols for c in late),
           f"missing {[c for c in late if c not in cols]} - header came from "
           f"{cols[1] if len(cols) > 1 else '?'}")
    built = rows.get("zz_built", {})
    record("statecols:and the measured row actually carries them",
           built.get("state") == "PRESENT" and (built.get("nt_expected_breadth") or ""),
           f"state={built.get('state')}, "
           f"nt_expected_breadth={built.get('nt_expected_breadth')!r}")
    record("statecols:the unbuilt row is UNKNOWN and simply blank there",
           rows.get("aa_unbuilt", {}).get("state") == "UNKNOWN",
           f"got {rows.get('aa_unbuilt', {}).get('state')}")



def tier_rank_local(tier: str) -> int:
    """Rank a tier for test comparisons, without importing the scorer's table."""
    order = ["no_reference", "screened_only", "attributed_to_neighbour",
             "indistinguishable_complex", "suppressed_control", "control_background",
             "marker_veto", "group_lead", "vaccine_strain", "supported",
             "confirmed_in_silico"]
    try:
        return order.index((tier or "").strip())
    except ValueError:
        return -1


def test_multi_reference_geometry(tmp: Path) -> None:
    """Breadth dilution is a property of having several references, not of rank.

    D42. `breadth` and `window_occupancy` are covered-over-SUMMED across every reference
    filed under a taxid, so a taxon carrying two assemblies has both halved by the one
    the sample does not contain. coverage_stats.py already emits
    `best_reference_breadth` / `best_reference_occupancy` on every row for exactly this,
    but score_evidence substituted them for supra-specific entries ONLY - so a species
    with two references was judged on a denominator including a genome it never had.

    The substitution is EXACT only when every assembly under the taxid is a single
    sequence, because `best_reference_*` is measured per BAM reference, i.e. per contig.
    With a draft assembly it is the best-covered contig, which must not stand in for a
    genome - so that case keeps the diluted numbers and says it may be under-called.
    """
    def refs(path: Path, rows: list[tuple[str, str, int]]) -> Path:
        path.write_text(
            "taxid\tname\tgroup\trole\taccession\tassembly_level\t"
            "genome_size\tn_contigs\tpath\n"
            + "".join(f"{tx}\tTestus {tx}\tbacteria\ttarget\t{acc}\tcomplete\t"
                      f"3000000\t{nctg}\t/x/{acc}.fna\n" for tx, acc, nctg in rows))
        return path

    # One row, two assemblies' worth of reference: the sample contains one of them, so
    # the summed denominator halves both statistics. The fixture leans on OCCUPANCY, not
    # breadth: a watchlist row's own `min_breadth` overrides `--threshold min_breadth`
    # (0.01 for this taxon), so a breadth floor cannot be set from here, whereas
    # `--min-occupancy` is a plain CLI floor. Diluted occupancy 0.476 fails a 0.60 floor;
    # 0.952 on the assembly actually covered clears it. So the substitution has to change
    # the TIER, not merely print a note.
    def row(sample: str) -> str:
        return cov_row(sample=sample, role="sample", platform="illumina", taxid="1639",
                       taxon="Listeria monocytogenes", group="bacteria",
                       reference_bp=5800000, unique_reads=80000,
                       reads_per_million=8000.0, breadth=0.4400,
                       mean_depth=40.0, depth_cv=0.40, window_bp=1000,
                       windows_expected=5800, windows_occupied=2761,
                       window_occupancy=0.4760, mean_identity=0.9900,
                       mean_alignment_length=99, best_reference_breadth=0.8800,
                       best_reference_occupancy=0.9520,
                       best_reference="taxid|1639|NC_003210.1",
                       best_neighbour_taxid="1642", best_neighbour_reads=20,
                       neighbour_ratio="4000.0", screen_assigned_reads=80000,
                       mapped_fraction_of_assigned="0.9900",
                       total_reads_screened=10000000, mapq_min=10)

    stats = tmp / "mref_stats.tsv"
    stats.write_text("file\tavg_len\tsample\trole\tplatform\n"
                     "/x/S.fq.gz\t101.0\tS\tsample\tillumina\n")

    def judge(tag, manifest_rows):
        cov = tmp / f"mref_{tag}.tsv"
        write_cov(cov, [row("S")])
        return score(tmp / f"d_{tag}", cov, markers=_empty_markers(tmp),
                     extra=["--refs-manifest", str(refs(tmp / f"m_{tag}.tsv", manifest_rows)),
                            "--read-stats", str(stats), "--min-occupancy", "0.60"])

    # (a) TWO single-sequence assemblies: the substitution is exact, so the row is
    #     judged on the assembly it is actually on and clears the floor.
    got = judge("two", [("1639", "GCF_000196035.1", 1), ("1639", "GCF_000008285.1", 1)])
    r = got.get(("S", "Listeria monocytogenes"), {})
    record("mref:a species with two single-sequence references is judged on the best one",
           "MULTI-REFERENCE TAXON:" in r.get("reasons", ""), r.get("reasons", "")[:150])
    record("mref:and the rescued row names the reference and both figures",
           "NC_003210.1" in r.get("reasons", "") and "0.8800" in r.get("reasons", ""),
           r.get("reasons", "")[:200])
    rescued_tier = r.get("tier")
    record("mref:the substitution changes the TIER, not just the prose",
           r.get("tier") not in ("screened_only", "no_reference"),
           f"got {r.get('tier')} - diluted occupancy 0.476 would fail the 0.60 floor")

    # (b) ONE assembly: nothing to dilute, so no note at all. This is the case that must
    #     not change - 167 of the pack's taxids are single-reference.
    got = judge("one", [("1639", "GCF_000196035.1", 1)])
    r = got.get(("S", "Listeria monocytogenes"), {})
    record("mref:a single-reference taxon gets no geometry note",
           "MULTI-REFERENCE TAXON" not in r.get("reasons", ""), r.get("reasons", "")[:150])

    # (c) TWO assemblies, one of them a draft: best_reference_* is a per-CONTIG figure,
    #     so substituting it would let one contig speak for a genome. Keep the diluted
    #     numbers and disclose it.
    got = judge("draft", [("1639", "GCF_000196035.1", 1), ("1639", "GCA_900000000.1", 58)])
    r = got.get(("S", "Listeria monocytogenes"), {})
    record("mref:a draft assembly blocks the substitution rather than risking a contig",
           "GEOMETRY NOT RE-MEASURED" in r.get("reasons", ""), r.get("reasons", "")[:150])
    record("mref:and the under-call is disclosed, not silent",
           "UNDER-called" in r.get("reasons", ""), r.get("reasons", "")[:200])
    record("mref:the blocked row really is held lower than the rescued one",
           tier_rank_local(r.get("tier")) < tier_rank_local(rescued_tier),
           f"blocked={r.get('tier')} vs rescued={rescued_tier}")

    # (d) No manifest: the question cannot be asked, so it must not be answered. This is
    #     the pre-D42 behaviour and every run without --refs-manifest depends on it.
    cov = tmp / "mref_nomanifest.tsv"
    write_cov(cov, [row("S")])
    got = score(tmp / "d_nomanifest", cov, markers=_empty_markers(tmp),
                extra=["--read-stats", str(stats), "--threshold", "min_breadth=0.50"])
    r = got.get(("S", "Listeria monocytogenes"), {})
    record("mref:without a refs manifest no geometry claim is made either way",
           "MULTI-REFERENCE TAXON" not in r.get("reasons", ""), r.get("reasons", "")[:150])



def test_refs_manifest_integrity(tmp: Path) -> None:
    """Two invariants over the shipped reference pack, both violated on 2026-10-21.

    D42. (a) A (taxid, accession) pair must appear ONCE. `--max-per-taxon N` deduped
    candidates by FILENAME, and a reused genome sits under the plain
    <taxid>_<name>.fna.gz while an extra is written as
    <taxid>_<name>.<accession>.fna.gz - so re-selecting the accession already on disk
    produced a different filename, passed the check, and was downloaded again. Running
    `--only-taxid 1639 --max-per-taxon 2` gave four taxa two manifest rows naming ONE
    accession and two byte-identical copies of one genome, while the summary reported
    "targets fully built: 1 / 1" and "genomes on disk: 8". N distinct assemblies were
    unobtainable for any taxon that already had one, i.e. every taxon in a built pack.

    (b) A competitive bundle must not be OLDER than the genomes it is built from. The
    FASTA and its minimap2 indexes are skipped whenever the file merely exists, so
    adding a reference left the manifest saying two and the aligner using one - a
    disagreement invisible from either side alone. Existence is not freshness.
    """
    refs = Path(os.environ.get("BTF_REFS_DIR",
                               "/media/8tb/bioThreatFinder/packs/human/refs"))
    man = refs / "manifest.tsv"
    if not man.is_file():
        skip("refs:every (taxid, accession) pair appears once",
             f"no reference pack at {refs} (set BTF_REFS_DIR)")
        return

    rows = [r for r in csv.DictReader(man.open(), delimiter="\t")]
    pairs: dict = {}
    for r in rows:
        key = ((r.get("taxid") or "").strip(), (r.get("accession") or "").strip())
        pairs[key] = pairs.get(key, 0) + 1
    dupes = sorted(k for k, n in pairs.items() if n > 1)
    record("refs:every (taxid, accession) pair appears once",
           not dupes, f"{len(rows)} rows; duplicated: {dupes[:6]}")

    # Distinct accessions, not just distinct rows: a taxon credited with N references
    # must hold N DIFFERENT assemblies, which is the thing the filename dedupe lost.
    per_taxid: dict = {}
    for r in rows:
        tx = (r.get("taxid") or "").strip()
        per_taxid.setdefault(tx, set()).add((r.get("accession") or "").strip())
    multi = {t: a for t, a in per_taxid.items() if len(a) > 1}
    counts: dict = {}
    for r in rows:
        tx = (r.get("taxid") or "").strip()
        counts[tx] = counts.get(tx, 0) + 1
    mismatch = sorted(t for t in per_taxid
                      if counts.get(t, 0) != len(per_taxid[t]))
    record("refs:a taxon's row count equals its distinct accession count",
           not mismatch,
           f"{len(multi)} multi-reference taxa; rows != distinct accessions for: "
           f"{mismatch[:6]}")

    # (b) freshness of the competitive bundle against its own target genomes.
    comp = refs / "competitive"
    if not comp.is_dir():
        skip("refs:no competitive bundle is older than its genomes",
             "no competitive/ directory")
        return
    stale = []
    for tx, _accs in sorted(per_taxid.items()):
        bundle = comp / f"{tx}.fna.gz"
        if not bundle.is_file():
            continue
        bt = bundle.stat().st_mtime
        for r in rows:
            if (r.get("taxid") or "").strip() != tx:
                continue
            gp = (r.get("path") or "").strip()
            if not gp:
                continue
            g = Path(gp)
            # A one-second slack: the bundle is written straight after the genome it
            # was just given, and identical mtimes are not staleness.
            if g.is_file() and g.stat().st_mtime > bt + 1:
                stale.append(f"{tx}:{g.name}")
                break
    record("refs:no competitive bundle is older than its genomes",
           not stale,
           f"bundle predates a genome it should contain: {stale[:6]}")



CASES = {
    "selftests": test_selftests,
    "tiering": test_tiering,
    "marker": test_marker_gating,
    "wright": test_wright_gates,
    "coverage": test_coverage_metrics,
    "screen": test_screen,
    "report": test_report_screen_coverage,
    "pack": test_pack_integrity,
    "depth": test_depth_gate,
    "eval": test_evaluator,
    "workflow": test_workflow_scripts,
    "platform": test_platform_carried,
    "absent": test_exhaustive_absent,
    "shared_assembly": test_shared_assembly_masking,
    "matrix_blank": test_matrix_blank_suppresses_host_artifact,
    "neighbour_named": test_neighbour_named,
    "contest": test_contest_attribution,
    "reachable": test_gate_reachability,
    "d17a": test_marker_block_required_to_confirm,
    "ctrlaudit": test_control_subtraction_audit,
    "crossbleed": test_cross_scenario_control_bleed,
    "protfloor": test_protein_tier_honours_catalogue,
    "markerexp": test_marker_expectation_cli,
    "tiergap": test_marker_tier_gap,
    "aacov": test_protein_coverage_per_sequence,
    "markerjoin": test_marker_catalogue_watchlist_join,
    "hopguard": test_index_hop_batch_guard,
    "defaults": test_scorer_defaults_match_shipped_config,
    "leadcap": test_supra_specific_entries_are_leads,
    "leadjoin": test_lead_join_seam,
    "ssrole": test_samplesheet_role_authority,
    "markergene": test_marker_gene_identity,
    "salm": test_salmonella_block,
    "statecols": test_marker_state_columns,
    "mref": test_multi_reference_geometry,
    "refs": test_refs_manifest_integrity,
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
