#!/usr/bin/env python3
"""check_pedestal_model.py - test `pedestal = recruitment_rpm x source_reads` against a run.

WHY THIS IS A SCRIPT AND NOT ARITHMETIC IN A COMMIT MESSAGE

The D12/D14 pedestal figure has been published four times and was wrong three times:
7.74% (records divided by pairs), 3.87% (right unit, measured at 101 bp for a 150 bp
batch), 5.51% (a 50k-read spot check quoted while a 500k run was producing a better
number), and a "within 6%" claim that was 6.31% on its own two inputs. Every one of those
was hand-derived from numbers copied between a terminal and a document.

So the comparison is computed here, from primary sources, and printed with the conditions
attached. Hand-copying a coefficient is the failure mode this file exists to retire.

WHAT IT READS
  recruitment matrix   scripts/measure_recruitment.py output (carries read_len and mapq)
  component JSON       <batch>/logs/<sample>.json - how many PAIRS of each organism were
                       simulated. Note `reads_requested`/`reads_actual` are PAIRS; the
                       pipeline's .nreads and total_reads_screened are RECORDS, i.e. 2x.
  coverage TSV         <run>/confirm/<sample>.<taxid>.coverage.tsv - the observed
                       best_neighbour_reads and unique_reads, in records.

WHAT IT REFUSES TO DO
  Invert the coefficient to estimate surviving HOST reads. Survivors of host depletion are
  not a random sample of host reads - they are the reads that failed to map to the host
  reference, which is the same class that cross-maps to a divergent eukaryotic genome. A
  coefficient measured on whole-genome simulated host reads is measured on a population
  the pipeline never sees, and inverting it produced an estimate larger than the sample's
  entire non-host remainder. The host arm is reported as an OBSERVED rate, never inverted.

    python3 scripts/check_pedestal_model.py \\
        --matrix $BTF_DB/benchmarks/recruitment_150bp_key.tsv \\
        --batch  $BTF_DB/simulated/lod_v1 \\
        --run    $BTF_DB/runs/lod_v2_kraken2 \\
        --source 562_Escherichia_coli.fna.gz --bundle 632 --neighbour 630
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys
from pathlib import Path


def die(msg: str) -> None:
    sys.exit(f"ERROR: {msg}")


def load_matrix(path: Path, source: str, bundle: str) -> dict[str, dict]:
    rows = [r for r in csv.DictReader(path.open(), delimiter="\t")
            if r["source"] == source and r["bundle_taxid"] == bundle]
    if not rows:
        die(f"no rows in {path} for source={source} bundle={bundle}")
    if "read_len" not in rows[0]:  # noqa: SIM102 - message is the point
        die(f"{path} has no read_len column, so the conditions it was measured under are "
            f"unknown. Regenerate with scripts/measure_recruitment.py - applying a "
            f"coefficient at the wrong read length is the error this check exists to stop.")
    return {r["ref_taxid"]: r for r in rows}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True, type=Path)
    ap.add_argument("--batch", required=True, type=Path, help="simulated batch dir")
    ap.add_argument("--run", required=True, type=Path, help="pipeline run dir")
    ap.add_argument("--source", required=True, help="matrix source filename")
    ap.add_argument("--organism", default="Escherichia coli",
                    help="component label in the batch JSON for --source")
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--neighbour", required=True)
    ap.add_argument("--bundle-taxon", help="taxon NAME for --bundle, e.g. 'Yersinia "
                    "pestis'. Without it a sample counts as spiked if it holds ANY spike, "
                    "which mislabels multi-scenario batches.")
    ap.add_argument("--target-arm", action="store_true",
                    help="Also compute the target arm, from the SAME samples the neighbour "
                         "arm used. The BAM set is derived, not supplied: an earlier "
                         "--attributed-bam let the operator choose, and the published "
                         "-0.6%% came from 4 of 12 available BAMs with no record of which. "
                         "Fixing the filter while leaving the subject set selectable moved "
                         "the laundering rather than removing it.")
    ap.add_argument("--organism-prefix", default="NC_002695.2,NC_002127.1,NC_002128.1",
                    help="comma-separated read-name prefixes identifying --organism reads")
    ap.add_argument("--samtools", default="samtools")
    a = ap.parse_args()

    refs = load_matrix(a.matrix, a.source, a.bundle)
    for t in (a.bundle, a.neighbour):
        if t not in refs:
            die(f"matrix has no row for ref {t} in bundle {a.bundle}")
    cond = refs[a.bundle]
    print(f"coefficients from {a.matrix.name}")
    print(f"  measured at read_len={cond['read_len']} mapq>={cond['mapq_min']} "
          f"insert_mean={cond.get('insert_mean', '?')} "
          f"over {int(cond['source_reads']):,} source reads")

    # Batch read length must match the coefficient's, or the comparison is meaningless.
    cmds = sorted(glob.glob(str(a.batch / "logs" / "*.cmds")))
    lens = set()
    for c in cmds:
        for line in open(c):
            if " -l " not in line or "art_illumina" not in line:
                continue          # only the ART line; a stray `ls -l` must not pollute it
            parts = line.split()
            i = parts.index("-l")
            if i + 1 < len(parts):     # the bounds check added for -m, missing here
                lens.add(parts[i + 1])
    if not lens:
        iss = any("iss " in open(c).read() or "iss generate" in open(c).read()
                  for c in cmds)
        die(f"no ART commands found under {a.batch / 'logs'}/*.cmds, so the batch's read "
            f"length is unknown."
            + (f"\n       This batch was simulated with InSilicoSeq, not ART. The "
               f"pedestal coefficients were measured under ART HS25 at a specified read "
               f"length AND insert size; iss has NO insert parameter and a model-fixed "
               f"read length, and insert alone moves coefficients 10-28%. A cliff measured "
               f"on an iss batch cannot be compared against an ART-calibrated prediction."
               f"\n       simulate_benchmark.py defaults --illumina-simulator to 'auto', "
               f"and 'auto' PREFERS iss when it is installed (it is). lod_v1 was built "
               f"with an explicit --illumina-simulator art. Rebuild with that flag."
               if iss else
               "\n       Previously this printed 'conditions match' and continued, which "
               "is worse than no guard: it asserts a condition it could not check."))
    if {str(cond["read_len"])} != lens:
        die(f"coefficient measured at {cond['read_len']} bp but the batch was simulated at "
            f"{sorted(lens)} bp. Recruitment is read-length dependent; this comparison "
            f"would be invalid.")
    # Insert size is checked too: the two 150 bp matrices on disk differ ONLY by insert
    # (300 vs 350) and that moves the target-arm coefficient by 10% and ref 633 by 28%.
    inserts = set()
    for c in cmds:
        for line in open(c):
            if " -m " not in line or "art_illumina" not in line:
                continue          # `samtools sort -m 4G` in the same log is not an insert
            parts = line.split()
            i = parts.index("-m")
            if i + 1 < len(parts):            # "-m" as the last token was an IndexError
                inserts.add(parts[i + 1])
    # Equality, not membership, and fatal when undeterminable - the read-length guard was
    # fixed this way and the insert guard was left with both holes: a batch with no -m at
    # all printed "conditions match", and mixed inserts passed because one of them matched.
    if not inserts:
        die(f"no -m (insert size) found in {a.batch / 'logs'}/*.cmds, so the batch's "
            f"fragment regime is unknown. Insert alone moves coefficients 10-28%.")
    if {str(cond.get("insert_mean", ""))} != inserts:
        die(f"coefficient measured at insert {cond.get('insert_mean')} but the batch used "
            f"{sorted(inserts)}. That alone moves coefficients by 10-28%.")
    print(f"  batch simulated at {sorted(lens)} bp, insert {sorted(inserts) or ['?']} "
          f"-> conditions match")

    rows = []
    no_organism, no_target, no_reads, no_coverage = 0, 0, 0, []
    for sample_json in sorted(glob.glob(str(a.batch / "logs" / "*.json"))):
        d = json.load(open(sample_json))
        sample = d["sample"]
        # Two different questions. The organism supplies the COEFFICIENT's denominator;
        # the bundle target is what stage 1 must nominate for a coverage TSV to exist. The
        # previous version tested the organism for both, so six negative controls - which
        # contain E. coli but never contained the target - were all labelled "stage 1 did
        # not nominate it".
        has_org = any(c.get("label") == a.organism for c in d["components"])
        # ANY spike is not THIS bundle's spike. On a multi-scenario batch (human_core.yaml
        # spikes Listeria, SARS-CoV-2, monkeypox, P. falciparum...) every such sample would
        # be labelled "stage 1 did not nominate it" for a bundle that was never in it - the
        # same mislabel fixed one level down for negative controls. --bundle-taxon narrows
        # it; without it the label is scoped to "a spike", and says so.
        spikes = [c.get("label", "") for c in d["components"] if c.get("kind") == "spike"]
        has_target = bool(spikes) if not a.bundle_taxon else a.bundle_taxon in spikes
        if not has_org:
            no_organism += 1
            continue
        if not has_target:
            no_target += 1
            continue
        cov = a.run / "confirm" / f"{sample}.{a.bundle}.coverage.tsv"
        if not cov.exists():
            no_coverage.append(sample)
            continue
        # reads_generated is what ART actually produced; reads_requested is the plan.
        # There is no reads_actual - an earlier version looked for it and so always
        # silently used the planned figure.
        # NOT an `or` chain. `reads_generated or reads_requested` silently falls back to
        # the PLAN whenever generated is 0 or absent - reinstating the 800,000 this script
        # exists to avoid, with no warning. ART under-produces routinely (B. thuringiensis
        # asked 100,000, got 98,241), so a low enough -f yields 0 and the fallback fires.
        comps = [c for c in d["components"] if c.get("label") == a.organism]
        if any("reads_generated" not in c for c in comps):
            die(f"{sample}: a {a.organism} component has no reads_generated. Falling back "
                f"to reads_requested would quote the plan as a measurement.")
        pairs = sum(c["reads_generated"] for c in comps)
        if not pairs:
            no_reads += 1
            continue
        src_reads = pairs * 2          # components are PAIRS; the pipeline counts RECORDS
        obs = next(csv.DictReader(cov.open(), delimiter="\t"))
        if obs.get("best_neighbour_taxid") != a.neighbour:
            # Comparing a coefficient for ref X against a count for ref Y is silent
            # nonsense, so it is refused rather than averaged in.
            die(f"{sample}: best_neighbour_taxid is {obs.get('best_neighbour_taxid')}, "
                f"not the --neighbour {a.neighbour} the coefficient was taken for.")
        rows.append((sample, src_reads,
                     float(refs[a.neighbour]["rpm_source_reads"]) * src_reads / 1e6,
                     int(obs["best_neighbour_reads"]),
                     float(refs[a.bundle]["rpm_source_reads"]) * src_reads / 1e6))
    # Separate reasons. The old count was (all JSONs - rows), so it swept in every sample
    # that never contained the organism at all and attached "stage 1 never nominated it" to
    # all of them.
    if not rows:
        die(f"no sample had both a component JSON and {a.bundle}.coverage.tsv")

    print(f"\n{'sample':<42} {'src reads':>10} {'pred nb':>8} {'obs nb':>7} {'err':>7}")
    errs = []
    for s, src, pred_nb, obs_nb, _ in rows:
        if obs_nb == 0:
            die(f"{s}: best_neighbour_reads is 0, so a relative error is undefined. This "
                f"is a real value - every negative 5811 row has it - not a parse failure.")
        e = pred_nb / obs_nb - 1
        errs.append(e)
        print(f"{s:<42} {src:>10,} {pred_nb:>8.0f} {obs_nb:>7,} {e:>+6.1%}")
    mean = sum(errs) / len(errs)
    print(f"\nNEIGHBOUR ARM: mean error {mean:+.1%} over {len(errs)} sample(s)")
    print(f"  skipped: {no_organism} with no {a.organism} component; "
          f"{no_target} with no spiked target at all (negative controls); "
          f"{no_reads} with the component but 0 generated reads; "
          f"{len(no_coverage)} spiked but with no {a.bundle} coverage TSV - stage 1 did not "
          f"nominate it: {', '.join(no_coverage) or 'none'}")
    total_json = len(glob.glob(str(a.batch / "logs" / "*.json")))
    accounted = len(rows) + no_organism + no_target + no_reads + len(no_coverage)
    if accounted != total_json:
        die(f"accounting does not close: {accounted} of {total_json} samples. A sample "
            f"vanishing from every counter is how 'skipped 3' came to be printed for 9.")

    # The target arm is reported too. Quoting only the arm that agrees is how a model gets
    # called calibrated when half of it is not.
    tgt_pred = sum(r[4] for r in rows) / len(rows)
    print(f"TARGET ARM:    predicted {tgt_pred:.0f} records on ref {a.bundle} per sample.")
    if a.target_arm:
        import subprocess
        prefixes = tuple(a.organism_prefix.split(","))
        obs_t, used = [], []
        bams = [a.run / "minimap2_competitive" / f"{s}.{a.bundle}.bam" for s, *_ in rows]
        missing = [b.name for b in bams if not b.exists()]
        if missing:
            die(f"target arm needs one BAM per neighbour-arm sample; missing: "
                f"{', '.join(missing)}")
        for bam in map(str, bams):
            refs = [ln.split("\t")[0] for ln in subprocess.run(
                [a.samtools, "idxstats", bam], capture_output=True, text=True,
                check=True).stdout.splitlines()
                if ln.startswith(f"taxid|{a.bundle}|")]
            if not refs:
                die(f"{bam}: no taxid|{a.bundle}| references")
            # -F 0x904 and -q are the SAME filters measure_recruitment.py uses. The
            # previous free-text values were -F 0x4, which includes supplementary
            # alignments the coefficient excludes.
            lines = subprocess.run(
                [a.samtools, "view", "-q", str(int(cond["mapq_min"])), "-F", "0x904", bam]
                + refs, capture_output=True, text=True, check=True).stdout.splitlines()
            obs_t.append(sum(1 for ln in lines if ln.startswith(prefixes)))
            used.append(Path(bam).name)
        m = sum(obs_t) / len(obs_t)
        if m == 0:
            die(f"no records matched --organism-prefix {a.organism_prefix!r} in any target "
                f"BAM. The default is O157:H7 Sakai's three accessions; a batch built on a "
                f"different E. coli assembly (e.g. K-12 NC_000913.3) needs its own. "
                f"Refusing to divide by zero - the neighbour arm got this guard first.")
        print(f"               observed {min(obs_t)}-{max(obs_t)} (mean {m:.1f}, "
              f"n={len(obs_t)}, same samples as the neighbour arm, counted here at "
              f"-F 0x904 MAPQ>={cond['mapq_min']}) -> {tgt_pred / m - 1:+.1%}")
    else:
        print("               NOT COMPARED: pass --target-arm. The coverage TSV does not "
              "separate which of the target's reads came from the background organism, so "
              "it is counted from the BAMs of the same samples under the same filter.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
