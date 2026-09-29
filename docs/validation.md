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
