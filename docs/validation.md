# Validation route

What it would take to make output from this pipeline defensible, and which
authorities apply. Assembled from a 2025–2026 literature and standards review; every
claim carries its source so it can be checked rather than trusted.

**Status: none of this has been done.** This document is the map, not the territory.

---

## The one document with enforceable numeric criteria

**New York State DOH, Wadsworth Center, CLEP: "Validation of Next Generation
Sequencing (NGS) Methods for Identification and/or Characterization of Infectious
Agents", March 2023.** Six pages, and the only document from any regulator that
explicitly addresses *metagenomic* assays with concrete acceptance criteria. Even
outside New York it is the right model, because it is the only one actually enforced.

What it requires, and where this repo stands:

| Requirement | Status here |
|---|---|
| LoD per **pathogen category** — DNA virus, RNA virus, gram-positive, gram-negative, mold, yeast, parasite | Scenarios cover all seven, tagged `nysdoh_category` |
| **6–8 ten-fold dilutions**, at least one to extinction | Decade series, 5–6 levels starting below any plausible LoD |
| **Separate LoD per specimen type** | Approximated by host fraction only. **Gap** — needs real matrices |
| **Co-infection, up to 4 agents**, mixed high/low | `coinfection_four_targets` |
| Accuracy: **30 clinical specimens** blinded vs comparator, or 10 near LoD if none | **Not started** — needs real specimens |
| Reproducibility: inter-assay ≥3 specimens × 3 days; intra-assay ≥3 × triplicate | **Not started** |
| **Single locked-down pipeline version** for validation | Partially — version and command are recorded; container digests are not pinned |
| Revalidate on **≥3 previously analysed runs** after any pipeline change | **Not implemented** |
| No-template control in **every** amplification step; process positive control; PhiX; internal/inhibition control when sequencing direct from primary specimen | Samplesheet has only `role=sample\|control`. **Gap** — see roadmap item B |
| SOP must describe database update procedures | **Not written** |

Two findings from it that are counter-intuitive and load-bearing:

- **"Specificity testing is not needed for non-targeted NGS applications intended to
  detect all DNA or RNA in a sample."** A regulator deliberately declining to impose
  a targeted-assay exclusivity framework on an untargeted one. Our specificity
  scenarios remain valuable — but they are about *near-neighbour misattribution*, not
  regulatory exclusivity, and should be described that way.
- **In-silico validation is explicitly endorsed:** "a pipeline's lower limit for
  accuracy can be assessed with **in-silico serial down-sampling** of sequencing
  reads. Simulation can also be used to test **intermixture of populations of species
  in the data at different levels of contamination**." That is a regulator describing
  this harness. `downsample_accuracy_floor` implements the first half.

---

## Normative stack

| Layer | Document | Note |
|---|---|---|
| Metagenomics-specific | **ISO 20397-3:2025** — Massively parallel sequencing, Part 3: metagenomics (2025-07-23, ISO/TC 276) | Covers database construction, bioinformatics, **and validation of analytical pipelines and databases**. The primary normative anchor. Also ISO/TS 24420:2023; whether Part 3 supersedes it is unverified |
| QMS frame | **ISO 15189:2022 cl. 7.3** (clinical) or ISO/IEC 17025:2017 | 7.3.3 forces full *validation*, not verification, for an in-house assay |
| Validation design | **CLSI MM17-Ed2** (multiplex nucleic acid assays) | Explicitly an "error-based" approach for assays where per-measurand validation is impractical — i.e. the mNGS problem. **There is no CLSI metagenomics document** (verified by enumerating the MM series) |
| LoD statistics | **CLSI EP17-A2E** + EP17IG | Assumes one measurand, so the workable adaptation is probit LoD **per representative organism class** — exactly what NYSDOH mandates |
| Dry lab | **Roy 2018** (AMP+CAP, J Mol Diagn 20:4) and **Duncavage 2023** (J Mol Diagn 25:3, with NIST/GIAB) | Duncavage is *the* citation for in-silico datasets substituting for physical samples in pipeline validation |
| Wet lab | **Schlaberg 2017** (ASM+CAP) | No formal successor from ASM/CAP/AMP/IDSA/ESCMID |
| Inspection hooks | CAP **MOL.36015** (wet bench), **MOL.36115** (bioinformatics), **MOL.36118** (LoD), **MOL.35880** (microbial nucleic acid contamination) | MOL.35880 was added in the 2022 edition and is the one explicitly metagenomics-driven item — the accreditation hook for contamination control |
| Pipeline documentation | **IEEE 2791-2020** BioCompute Object | FDA-collaborated machine-readable format |
| Reporting | **STROBE-metagenomics** (Lancet Infect Dis 2020;20:e251) | |
| SOP scaffolding | **CDC/APHL NGS Quality Initiative** — 113 free QMS tools | Names "agnostic (metagenomic) sequencing" as a gap where tools do not yet exist |

### Animal health
**WOAH Terrestrial Manual Ch. 1.1.7** (HTS/bioinformatics; revised May 2025, 92nd
General Session Res. 26) plus **Ch. 1.1.6** and 2.2.1–2.2.8. Terrestrial only — no HTS
chapter in the Aquatic Manual. Three provisions matter directly:

- SOPs for the whole chain must exist **before validation can start**, and **Stage 1
  data must confirm analytical sensitivity, specificity and repeatability.**
- **"Confirmatory adjunct or secondary diagnostic assays need to be validated only
  for their analytical performance … and not to the full diagnostic extent (DSe/DSp,
  stage 2)."** This formally sanctions a screen-then-confirm architecture.
- **Where depth of coverage falls below the validated minimum, results should be
  confirmed by an independent method (e.g. Sanger).** The clearest normative
  confirmation requirement from any body.

The Aquatic Manual's per-method **"LV = level of validation"** grading is a useful
precedent for *grading* an untargeted method rather than passing or failing it.

### Plant health
**EPPO PM 7/151 (1)** (HTS in plant health diagnostics, 2022, + 2024 addendum),
**PM 7/98 (6) 2025** (accreditation; validation lives here), **PM 7/122 (3) 2025**
(interlaboratory comparisons). EPPO has the most mature untargeted-method guidance of
any body. Methodological literature: Massart 2022 (Peer Community J 2:e62), Lebas
2022 (EPPO Bull 52:394), and **Tamisier 2021** (Peer Community J 1:e53) on
semi-artificial datasets for validating plant-virus detection pipelines — the plant
analogue of this harness. IPPC **ISPM 27** has only generic validation language and
nothing on HTS.

---

## Reference materials

| Material | What it gives | Caveat |
|---|---|---|
| **NIST RM 8376** | 19 bacterial gDNAs + human gDNA, **ddPCR copy-number traceable** | The only metrologically anchored pathogen-DNA panel for mNGS. DNA-only, so it does not exercise extraction or lysis bias. **Apparently discontinued** |
| **NIST RM 8048** | Human faecal material, omnivore + vegetarian (issued 2025-03-07) | The only real-matrix NIST RM in this space |
| **ZymoBIOMICS D6310/D6311** | Community Standard II, **log-distributed 10²–10⁸ cells** | Purpose-built for LoD determination — the item you actually want |
| ZymoBIOMICS D6300 / D6305 / D6320-D6321 | Whole cells / DNA / spike-in controls | Match the `spike_in` enums UKHSA mSCAPE uses |
| **ATCC MSA-4000** | gDNA from 11 clinically common/AMR strains, ddPCR-quantified | **Independent work reports label-vs-actual species discrepancies for some ATCC Microbiome Standards — verify before using as a truth set** |

## External quality assessment

**QCMD "2026 Viral Metagenomics NGS EQA Pilot Study" (NGSmeta_26)** — 5 panel
members, synthetic CSF + human cell lines mimicking enterovirus/HSV/influenza, scored
at family/genus/species/subtype. Ran identically in 2024 and 2025.

It is the **only** untargeted-metagenomics EQA that exists, it is still labelled a
**pilot**, and the ISO 17043 accreditation statement that appears on QCMD's mature
schemes is **absent** from it. No GenQA, UK NEQAS, ECDC, WHO or CAP equivalent.
CAMI is benchmarking, not EQA — no panel material, no per-lab certificate.

**Document the PT gap as a QMS non-conformance/risk rather than ignoring it.** CAP's
in-silico NGS bioinformatics surveys (NGSB\*) are the nearest structural analogue and
have a four-year performance record (Furtado 2023, Arch Pathol Lab Med 147:2).

---

## US regulatory position, as of September 2026

- The **2016 FDA draft guidance on infectious-disease NGS was withdrawn on
  2023-04-05** and never finalised. There is no successor. The in-silico and
  FDA-ARGOS validation concessions it floated **have no current regulatory standing.**
- The **LDT final rule was vacated nationwide on 2025-03-31** (*ACLA v. FDA*, 776 F.
  Supp. 3d 554, E.D. Tex.) and FDA formally rescinded it at 90 Fed. Reg. 45134
  (2025-09-19), reverting 21 CFR 809.3(a).
- So clinical mNGS remains an **LDT under CLIA (42 CFR 493)** plus a CLIA-deemed
  accreditor (CAP) plus NY CLEP for New York specimens. **No FDA premarket pathway is
  being applied, and there is no FDA analytical-validation expectations document in
  force.** Verified against openFDA: no 510(k), PMA or classification record matches
  "metagenomic".
- **FDA-ARGOS** (BioProject PRJNA231221) inclusion criteria are worth mirroring in our
  own reference-set policy: organism identified by an independent reference method
  before sequencing; **two orthogonal sequencing technologies**; de novo assembly with
  **≥95% of positions at ≥20× depth**.

---

## When orthogonal confirmation is actually required

In descending order of authority:

1. **WOAH 1.1.7** — confirm by an independent method where **coverage depth falls
   below the validated minimum**. Also recognises a reduced validation tier for
   confirmatory assays, i.e. the framework assumes sequencing can itself *be* the
   confirmatory test.
2. **NYSDOH CLEP** — requires a documented **discrepancy-resolution procedure** and
   flagging of incidental findings unsupported by peer-reviewed evidence. Not blanket
   per-call confirmation.
3. **CAP** — confirmatory testing is among the required elements; a defined
   confirmation policy is expected.
4. **Professional guidance** (Schlaberg 2017, ENNGS Part II, Kan 2024) — confirmation
   framed around validation and around unexpected or actionable findings.

**Defensible formulation:** orthogonal confirmation is required (a) during validation
as an accuracy comparator, (b) for any detection below validated coverage/read
thresholds, and (c) for unexpected, novel or clinically actionable findings. **No body
requires per-call confirmation of every organism in every report** — so the blanket
"requires orthogonal laboratory confirmation" language this pipeline stamps on every
finding is *more* conservative than any standard demands. That is a defensible
posture for a biothreat tool, but it should be a deliberate choice, not an accident.

---

## Honest gaps

- **No specimens, no wet lab, no comparator data.** Everything above stage
  "in-silico" is unstarted.
- **Specimen-type LoD** is approximated by host fraction. Real matrices are required.
- **Container digests are not pinned**, so "single locked-down version" is not yet
  satisfiable.
- **No revalidation-on-change mechanism** (NYSDOH: ≥3 previously analysed runs).
- **No published method exists** for vaccine-vs-wild-type calling from metagenomes,
  nor for *Y. pestis* vs *Y. pseudotuberculosis* or *F. tularensis* subspecies
  discrimination from metagenomic data. Our marker tier addresses real gaps, which
  also means it is unvalidated and must not be presented as established practice.
- **STRATUM** (*J AOAC Int* 2026, doi:10.1093/jaoacint/qsag085) is a multi-agency,
  consequence-tiered reporting framework built for exactly this use case, including
  national-security and synthetic-biology detection. **Not yet read.** It should be,
  before the reporting layer is considered settled.

## Measured detection limit — 2026-10-03

First measurement against ground truth, `bin/evaluate_benchmark.py` over
`benchmarks/scenarios/human_core.yaml` and the `lod_v1` manifest (21 samples, 3
replicates per level).

**LoD = 1000 RPM** for *Y. pestis* on Illumina, and it is a cliff rather than a gradient:
0.000 detection at 1, 10 and 100 RPM, 1.000 detection AND 1.000 confirm at 1000 and
10000. `lod_half_rpm` equals `lod_all_rpm`, so there is no partial-detection band to
interpolate. Flagged `coarse - fewer than 5 replicates`.

**This is the first empirical read on decision D2.** sylph-only (`runs/lod_v1`) and
Kraken2 (`runs/lod_v2_kraken2`) produce the SAME LoD and the same 6
`confirmed_in_silico` findings at PPV 1.000. Kraken2's distinct-minimizer evidence did
not lower the detection limit on this batch. It is not a null result for D2 overall -
this batch has one target, one platform and no divergent or novel agent - but it does
mean the primary arm is not buying sensitivity here.

**What Kraken2 did change, and it needs calibration work:**

| tier | sylph-only | Kraken2 |
|---|---|---|
| `confirmed_in_silico` | 6, PPV 1.000 | 6, PPV 1.000 |
| `attributed_to_neighbour` | 3, PPV 1.000 | **68, PPV 0.500** |
| `screened_only` | — | 47, PPV 0.000 |
| `supported` | — | 18 (all taxa absent from benchmark truth, so unassessable) |

`attributed_to_neighbour` counts as NEGATIVE in the evaluator's mapping, so those 6
truth-present rows are scored as missed detections. One scenario expectation fails in
BOTH runs for the same reason:

    FAIL lod_gram_negative/Yersinia pestis min_tier_at_rpm[100]:
         expected >= supported, observed attributed_to_neighbour

At 100 RPM *Y. pestis* carries ~2,300 unique reads at 0.062 breadth and is reassigned to
a neighbour. Either the scenario's expectation is optimistic or `min_neighbour_ratio` is
miscalibrated; the measurement cannot tell these apart, and resolving it is step 4.

**Specificity, as measured:** 0.000 false-positive and 0.000 false-confirm across all 6
negative-control samples - **for the taxa the scenario lists**. That qualifier turned out to
matter: see D14. `negative_host_only` names eleven taxa in `truth.absent`, and *Toxoplasma
gondii* is not one of them, so the 18 `supported` *T. gondii* calls in this batch - every one
driven by 100% human reads - were scored `unlisted` and excluded from the specificity table
entirely. The 0.000 is real for *B. anthracis* and *F. tularensis*; it is not a pack-wide
specificity figure, and this file should not have implied otherwise.

## Kraken2 memory mapping — verified 2026-10-03

Two samples, `params.kraken2_memory_mapping = true`, `maxForks = 2`:

| | no mmap | mmap |
|---|---|---|
| KRAKEN2 realtime | 30 s - 3 m 24 s | **11.0 s / 11.1 s** |
| peak RSS | 115.5-115.7 GB | 115.3 GB |
| concurrent tasks | 1 | **2, overlapping** |

Both tasks were submitted 13 s apart and ran at once with peak RSS FLAT, which is the
point: they shared one file-backed copy of the 123.5 GB hash in page cache. Without mmap
the same pair would have demanded ~247 GB of unreclaimable memory.

## D12 mechanism confirmed by read attribution — 2026-10-03

The simulated reads carry their source accession in the read name, so the competitive
alignment BAMs already on disk answer the D12 question without any new sequencing. For
`runs/lod_v2_kraken2`, taxid 632 bundle, MAPQ>=10:

**Neighbour reference (taxid 630, *Y. enterocolitica*) — 6 samples spanning 4 spike levels:**
2,708–2,822 reads assigned, of which **97.4–98.2% originate from *E. coli*** O157:H7 Sakai
(the 2% standing background) and **0 originate from the spiked *Y. pestis*, in every
sample**. The rest is *H. haemolyticus* and *N. lactamica*, also background.

**Target reference (taxid 632) — the numerator is affected too:**

| level | target reads | true *Y. pestis* | *E. coli* | % false |
|---|---|---|---|---|
| 10 rpm | 427 | 225 | 202 | 47.3% |
| 100 rpm | 2,336 | 2,114 | 222 | 9.5% |
| 1000 rpm | 22,320 | 22,128 | 191 | 0.9% |
| 10000 rpm | 221,439 | 221,209 | 230 | 0.1% |

Both pedestals are flat in absolute terms (~2,700 on the neighbour, ~210 on the target) and
their ~13x asymmetry is a property of the references, not of the specimen. Consequence:
`neighbour_ratio ~= (true_target + 210) / 2,700`, so the 1.5 threshold needs ~3,840 true
target reads whether or not a neighbour is present.

**Prediction for the finer grid:** 2,114 true reads at 100 rpm scales linearly, so 3,840 is
reached at **~180 rpm**. The cliff should sit there, not at 1000 - the 1000 figure is an
artefact of a series that steps 100 -> 1000. Finding it elsewhere falsifies the arithmetic.

At 10 rpm nearly half the "target" reads are *E. coli*, so breadth, `depth_cv` and
`mean_identity` are all affected at low abundance, not just the ratio.

Masking cannot currently help: `mask_pack_refs.py`'s `crossgroup` source compares only
across distant groups by design, so bacteria-vs-bacteria recruitment is 0 bp masked for both
references. That makes within-family masking an untried lever.
