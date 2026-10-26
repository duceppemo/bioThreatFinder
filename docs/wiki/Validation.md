# Validation

**Read this before drawing a conclusion from a run.**

This pipeline is research software under active development, **not a validated assay**. The
gap is specific, and it is documented rather than glossed.

## What is not established

- **No threshold in `packs/*/pack.yaml` is a measured operating point.** They are reasoned
  defaults. Three are explicitly marked PROVISIONAL in the code, and
  `scripts/check_threshold_coherence.py` reports **75 of 89** taxa carrying a
  `min_unique_reads` that can never bind — up to 243× too high for *Schistosoma mansoni*.
- **Sensitivity and specificity are not characterised for any intended use.** What has been
  measured is below; it is a long way from a performance claim.
- **The marker tier addresses real gaps in the published literature, which is also why it
  is unvalidated** and must not be presented as established practice. There is no published
  framework for vaccine-versus-wild-type calling from metagenomic data, none for
  *Y. pestis* / *Y. pseudotuberculosis*, and none for *F. tularensis* subspecies
  discrimination from metagenomes.
- **The watchlist is not an authoritative legal source** for any jurisdiction.
- **There is no CLSI standard for untargeted metagenomic sequencing**, and no accredited
  mNGS proficiency scheme was found by enumeration of the MM series.
  [`docs/validation.md`](https://github.com/duceppemo/bioThreatFinder/blob/main/docs/validation.md)
  sets out the route to a defensible result and where the normative gaps are.

## What has been measured

- **Real data, Illumina and ONT.** ZymoBIOMICS D6300 (EVEN and LOG) plus all ten member
  organisms as single-organism libraries, 54 GB, md5-verified. It refuted a neighbour gate
  that had been reasoned rather than measured (now defaulted off) and re-scored the tier
  ladder, cutting real-data false positives **12 → 2** with both simulated benchmarks
  unchanged.
- **A limit of detection**, on one organism in one simulated matrix: 1000 RPM for
  *Y. pestis* on Illumina, with the cliff located to 176.6 rpm on a fine grid.
- **`confirmed_in_silico` has fired on real data**, backed by curated markers — *S.
  enterica* in a public library, markers at breadth 1.000 and identity 0.985–0.989, exactly
  one row of 82 moving.
- **The first ONT run end to end**: 15h 36m, 425 CPU-hours. Specificity clean — 84 absent
  taxa, four findings, all held at `attributed_to_neighbour`, **zero** false
  confirmations. Sensitivity one of four.
- **A measured negative-control effect.** Declaring one library a control moved `supported`
  PPV **0.600 → 0.750** and false positives 2 → 1, with no per-scenario detection rate
  changing.

## The gaps that matter most

**ONT has never been tested on the chemistry its floor is calibrated for.** Every ONT call
in this project is capped by the D10 identity floor and nothing else: *L. monocytogenes*
with breadth 0.95, 136,951 rpm, neighbour_ratio 9.7 and **both** confirmatory markers
PRESENT failed on `identity 0.8579 >= 0.9000` alone. At 0.90 no ONT sample can reach
`confirmed_in_silico` at all. **Do not lower the floor on this evidence:** the data is
R9.4.1, and 0.90 was reasoned for R10.4.1 at ~97–98%. The floor was not falsified; this
benchmark cannot test it.

**19 of 90 watchlist taxa cannot be reported by the Kraken2 arm at all** — no node in its
pruned taxonomy. All 19 are reachable through sylph, so one arm is blind rather than the
stage, but for that group sylph is not a redundant second opinion: it is the only arm that
can see them.

**31 of 90 taxa have no marker rule**, so `confirmed_in_silico` is unreachable for them by
construction. That cap is correct — the binding constraint on the top tier is curation, not
a threshold.

**Four near-neighbour genomes do not exist anywhere** (*Giardia ardeae*, *Entamoeba coli*,
Sandfly fever Naples phlebovirus, rat hepatitis E virus), recorded with empty paths so the
weakened discrimination is visible rather than silent. And *F. tularensis* versus
*F. novicida* cannot currently be discriminated: both resolved to the same assembly.

**Still missing entirely:** proficiency panels, characterised clinical specimens, a declared
no-template control, and any ONT data on R10.4.1.

## A negative is weaker than it looks

Almost every bug found in this project shares one shape: *an empty or degraded result
indistinguishable from a clean negative.* Several guards exist specifically to make that
visible — the stage-1 coverage table, the database freshness checks, the three-valued
marker logic — and none of them can turn absence of evidence into evidence of absence.

When reading a clean findings table, check the stage-1 coverage table and the unclassified
fraction too ([Outputs](Outputs)).

## The regression suite

```bash
python3 tests/test_pipeline.py          # 512 checks, offline, no databases, ~30 s
python3 tests/test_pipeline.py --list   # the groups
```

It asserts **designed behaviour, not accuracy** — a failure means the logic changed; it does
not mean the pipeline got less accurate, because nothing in it measures accuracy. On a bare
clone 446 pass and 16 skip, naming what to set. CI runs it on two Python versions, compiles
`main.nf`, and builds and smoke-tests the container images.

New checks are **mutation-tested**: break the thing the check is for, and watch it fail.
Three checks in this repository passed a mutation they should have caught and were
rewritten. Equally, four guards were caught crying wolf on a *healthy* installation and
fixed before shipping — a guard that fires on a correct system gets switched off, which is
worse than not having it.

## Benchmarking

32 scenarios, 29 machine-checkable expectations, weighted deliberately towards failure
modes rather than easy positives: near neighbours with the target genuinely absent, target
chromosomes stripped of their virulence plasmids, vaccine strains, index hopping, a
99.9%-host low-biomass case, and a plumbing positive control whose only job is to fail
loudly on misconfiguration. Limit-of-detection design follows NY State DOH CLEP (March
2023) — all seven pathogen categories, decade dilutions with one to extinction, and
four-agent co-infection.

`no_reference` observations are excluded from every rate, so adding watchlist taxa without
building their references cannot fake an improvement in specificity.

Simulation flatters every tool and establishes relative performance only, not clinical
performance.
[`docs/benchmarking.md`](https://github.com/duceppemo/bioThreatFinder/blob/main/docs/benchmarking.md)
states what it cannot tell you, and the two limits benchmarking never removes: software
cannot beat its input, and homology-based detection cannot see a wholly novel agent.

## Where the detail lives

- [`docs/verification.md`](https://github.com/duceppemo/bioThreatFinder/blob/main/docs/verification.md)
  — what has been checked, how to re-run it, and an explicit list of what has **not**
- [`docs/validation.md`](https://github.com/duceppemo/bioThreatFinder/blob/main/docs/validation.md)
  — normative standards, reference materials, EQA, regulatory position
- [`docs/benchmarking.md`](https://github.com/duceppemo/bioThreatFinder/blob/main/docs/benchmarking.md)
  — benchmark design and the tool bake-off
- [`docs/roadmap.md`](https://github.com/duceppemo/bioThreatFinder/blob/main/docs/roadmap.md)
  — what is not built, ranked by published effect size
- [`docs/decisions.md`](https://github.com/duceppemo/bioThreatFinder/blob/main/docs/decisions.md)
  — every choice, its evidence tag, and what would reverse it. Including the ones that were
  reversed
