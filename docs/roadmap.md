# Roadmap

## Blocking a real run
1. ~~Upgrade Nextflow~~ DONE - Nextflow self-updated to 26.04.6, which satisfies
   the manifest requirement. Note this forced a strict-syntax migration:
   `check_max()` config functions -> `resourceLimits`, no top-level `def` in
   config, no top-level `workflow.onComplete`, and process directives must precede
   `input:`. All migrated.
2. ~~Create the tool env~~ DONE - and it had to become **three** envs, not one.
   Installing kraken2 alongside mmseqs2/hmmer/mafft made the solver resolve
   kraken2 2.0.8-beta (2019), which has no `--report-minimizer-data`. `btf-refs` holds
   datasets/minimap2/seqkit/samtools/mmseqs2/hmmer/mafft; `btf-kraken2` holds
   kraken2 2.17.2 and bracken; `btf-sylph` holds sylph 1.0.0, matching the pin in
   `modules/local/detect.nf` because a sketch must be readable by the version that
   queries it. See docs/STATUS.md.
3. ~~Build databases~~ DONE - taxonomy snapshot, pack refs (90 targets, 379 rows),
   masking, markers, host genome, sylph pack sketch (rebuilt 2026-10-21 over all 375
   references), Kraken2 PlusPF and Metabuli are all built on NVMe. KrakenUniq is
   deliberately not built (D2). Freshness is checked for the competitive bundles (D45) and
   for the sylph sketch (D46 — recorded at build time in `<db>.genomes.txt` and compared
   against the pack at launch, plus `scripts/build_sylph_db.sh --check`); the Kraken2,
   Metabuli and marker databases still notice nothing when the references underneath them
   change.
4. **Validate on known-positive data** before trusting any output. *Partly done:*
   ZymoBIOMICS D6300 (EVEN and LOG, Illumina and ONT) plus all ten member organisms have
   run end to end, which refuted D12's neighbour gate and re-scored the tier ladder
   (D18, D35, D38). Still missing: proficiency panels, characterised clinical specimens,
   a declared no-template control, and ONT data on R10.4.1 chemistry. **No threshold in
   `packs/*/pack.yaml` is a measured operating point yet.**

## Next feature, in priority order

### DONE - Marker-gene / plasmid module
Implemented in v0.2.0: 166 markers, 59 rule sets, three-valued rule engine,
three evidence tiers (dna/protein/hmm). Every high_consequence *species* has a rule;
the only high_consequence row without one is the `Filoviridae` group entry, which is
capped at `group_lead` and so has nothing for a rule to gate.
See README for the design and validation. Remaining marker work:
  * 31 of 90 taxa still have no rule — 28 notifiable, 2 emerging, and Filoviridae.
    *Orientia tsutsugamushi* is the last notifiable bacterium among them
  * curate_required markers: B. mallei absent regions, vaccinia-specific loci,
    F. novicida discriminating loci, C. auris clades
  * absence-region logic as a first-class marker class. B. mallei needs it, and
    so does `bruc_S19_eri` - D30 established that the S19 marker CANNOT be fixed
    by retrieval, because the eri operon's presence is the wild-type state and
    supplying it would make `vaccine_if` fire on wild-type B. abortus
  * 4 retrieval gaps left after D30, each with a measured diagnosis in the
    catalogue note: `henipa_cedar` (gene deposited as bare `N` under
    *Henipavirus cedarense*; needs an accession), `orthopox_ATI` and `mpox_C3L`
    (need the ORF/OPG designation as an alias - the records are whole genomes
    whose CDSs use poxvirus ORF names), `lmono_LIPI4` (CC4-specific island,
    absent from EGD-e, so the lmo* tags are the wrong reference)
  * 9 `motif_todo` markers are now the dominant blocker: they alone make 9 of
    10 `vaccine_if` rules inoperative

### 1. Benchmark harness (now the highest value - nothing is measured)
Without ground truth, "accurate" is not assertable and every threshold is a guess.
  * Simulated reads: ART / InSilicoSeq (Illumina), NanoSim / Badread (ONT)
  * Deliberate near-neighbour spike-ins: B. cereus for anthracis, B. thailandensis
    for pseudomallei, E. dispar for histolytica, vaccinia for variola, Reston for
    Zaire, B. holmesii for pertussis
  * Vaccine-strain spike-ins: BCG, vaccinia, Sabin, 17D, RB51, TC-83
  * Mock communities (ZymoBIOMICS) and a spiked dilution series for limit of
    detection per taxon
  * Report sensitivity, specificity and LOD per taxon; commit the results

### 2. Empirical null and calibrated FDR
Replace fixed cutoffs with per-taxon background distributions from a panel of
negative and unrelated metagenomes. Report a calibrated statistic alongside the
interpretable gate breakdown.

### DONE - Stage-1 classifier upgrade (v0.3.0)
Kraken2 with --report-minimizer-data is the primary arm; KrakenUniq is opt-in and off
by default (this was reversed after reading Kraken2's source and Oskolkov 2026 - see
docs/decisions.md D2). Metabuli supplies a DNA+amino-acid arm for divergent agents and is
likewise opt-in and off by default since 2026-10-05 (D15 - it changed no finding on any
scenario measured, including the divergent_* pair built to favour it). Bracken abundance,
sylph ANI. Any-classifier OR; concordance and k-mer quality recorded and
surfaced in the report. Project DBs moved to NVMe at /media/8tb/bioThreatFinder.
Remaining stage-1 work:
  * Centrifuger as an ONT-specific arm (low memory, FM-index) - evaluate whether
    it adds anything over KrakenUniq on long reads
  * Feed distinct-k-mer evidence into the empirical null rather than only flagging

### 3. Accuracy upgrades to stage 2
  * Stage 2: multi-reference / pangenome bundles per taxon instead of one
    reference, with skani for ANI-based assignment
  * Assembly tier: metaSPAdes / Flye on candidate reads, then confirm on contigs

### 4. Containerisation
Pinned digests per process, Apptainer on HPC. Required before any result is
defensible or publishable.

### 5. Old priority - vaccine / attenuated-strain discrimination
For several watchlist organisms the chromosome is NOT the actionable signal, so the
current pipeline cannot make a defensible call:

| Organism | Actionable locus | Why the chromosome is not enough |
|---|---|---|
| *B. anthracis* | pXO1 / pXO2 | Chromosome is near-identical to the B. cereus group |
| *C. botulinum* | BoNT cluster | Species is polyphyletic; cluster also in *C. butyricum*/*baratii* |
| STEC | stx1 / stx2 / eae | *E. coli* presence is uninformative |
| *V. cholerae* | ctxAB / tcpA | Environmental non-toxigenic strains are common |
| *C. diphtheriae* | tox | Phage-borne, not chromosomal |
| *B. pertussis* | discriminating loci | IS481 is multi-copy and shared with *B. holmesii* |
| *Y. pestis* | pPCP1 / pMT1 | >97% chromosomal ANI to *Y. pseudotuberculosis* |
| *L. monocytogenes* | LIPI-1 | Discriminates from *L. innocua* |
| *A. fumigatus* | cyp51A alleles | Resistance, not presence, is actionable |

Design: add a `markers/` directory per pack keyed by watchlist taxid, align candidate
reads against the marker set, and make `confirmed_in_silico` conditional on marker
coverage where the watchlist demands it. This needs a new `requires_marker` column.

Partially DONE: vaccine_strain markers and a `vaccine_if` rule exist for BCG,
RB51/S19, Sabin, 17D, Candid#1, TC-83, MP-12/Clone 13, measles Edmonston and
rabies ERA/SAD. Still needed: curated vaccinia-specific loci (the dominant
variola false positive), and explicit competitive references for each vaccine
strain rather than marker motifs alone.

### 3. Read extraction for candidates
An analyst cannot currently eyeball alignments. Emit per-finding BAM slices plus a
coverage plot.

### 4. Animal pack (WOAH listed diseases)
Multi-host depletion (bovine, porcine, ovine, caprine, equine, avian, salmonid).
Reconcile the watchlist against the current WOAH list.

### 5. Plant pack (EPPO A1/A2)
Heavier on fungi/oomycetes. Crop host genomes. Note `/media/36tb/db/oomycetes` and
`stagonospora` already exist on this box.

### 6. Influenza A segment-aware subtyping
HA/NA subtype is the actionable unit; a species-level call is not enough. Needs a
dedicated segment-aware caller, flagging H5/H7 explicitly.

### 7. AMR characterisation as stage 3
CARD, ResFinder and DeepARG databases are already present on this box.

## From the 2025-2026 literature review (not yet implemented)

Ranked by published effect size. Sources in each line.

### 0. CORRECTION to earlier advice in this file and in conversation
I previously told the user target enrichment was the single biggest sensitivity lever,
"orders of magnitude" ahead of software. **For a mission that must catch DIVERGENT
agents that is wrong, and measurably so.** Yamaguchi et al. 2025 measured it:

  * Capture needs **>65% nucleotide identity over >=1 kb contiguous** to enrich at all.
  * Off-panel agents get WORSE: simian haemorrhagic fever virus read abundance fell
    **25-fold** after capture; Nodamura virus segment 1 fell 33-fold, segment 2 50-fold.
  * A divergent bunyavirus that assembled readily from shotgun reads was "orders of
    magnitude lower" in the capture library.
  * Panels also age: probes designed on 2017-2018 influenza gave 40-86% coverage for
    H5N1 clade 2.3.4.4b vs >98.5% for a 2024-designed panel.

So enrichment is a large lever for ON-PANEL targets in low-biomass matrices and an
actively NEGATIVE one for the novel agent this tool exists to catch. The published
answer if both are wanted is **teamNGS**: pool enriched and unenriched libraries
90:10 or 80:20 in one run, which beat mNGS alone (p=0.0001 viral) and capture alone
(p=0.02). **Always keep an unenriched arm.**

Adaptive sampling is likewise not the lever I implied: 1.2-5x realised enrichment,
50-80% yield cost, 400 bp fragment floor. Treat it as free polishing, not strategy.

### 0b. THE novel-agent failure mode - act on this first
Perseus (Bioinformatics 2026, doi:10.1093/bioinformatics/btag687) measured what
Kraken2 does with sequence whose family is absent from the database: **the
false-assignment rate approaches 50% on long reads.** Divergent reads are not left
unclassified, they are CONFIDENTLY ASSIGNED TO THE WRONG KNOWN TAXON. Species
assignments were frequently supported by a median of only five k-mers.

For surveillance this is the worst possible shape of error: the novel agent is not
missed, it is mislabelled as something familiar and dismissed as background. No
amount of stage-2 confirmation on the WRONG taxon recovers it.

Partially addressed: the unclassified fraction is now a first-class per-sample output
with a HIGH_UNCLASSIFIED flag, because it is the only signal that survives this
failure mode. Still to do:
  * a Perseus-style abstention/rank-back-off layer, so the honest answer "genus
    Bacillus, insufficient evidence for anthracis" can be emitted
  * triage of the unclassified fraction (assemble it, protein-search it) rather than
    only counting it
  * ~~keep Metabuli despite its poor independent precision - its role here is RECALL on
    divergent sequence, and stage 2 supplies the precision.~~ **Settled 2026-10-05 by
    measurement, and not the way this bullet assumed.** The recall argument was never
    tested when it was written; it has been now, on the two `divergent_*` scenarios built
    for exactly this question, and the arm recovered nothing the nucleotide arm missed -
    both assign all 3,992 spiked reads to the same taxon, read for read. Across every
    scenario run it changed zero findings at ~22x Kraken2's CPU. Kept but defaulted OFF,
    because its one untested case - an agent absent from the index - cannot be simulated
    from a database. See D15.

### A. Mask the reference database before tuning any threshold
Now confirmed by two independent reviews. VirMask-style removal of recurrent
contaminants (adapters, PCR controls, bacterial rRNA/gDNA, CMV enhancers, phage) cut
false detections **40-94% with 100% true-positive preservation** (J Clin Microbiol
2026) - described in the second review as a better specificity trade than anything
available downstream, including assembly; Conterminator-style cross-contamination screening is described as
table stakes. Independently, Wright et al. 2026 traced an *H. pylori* false positive
to **phage sequence inside the reference genome itself**.

**DONE 2026-09/10.** `scripts/mask_pack_refs.py` is implemented, run and shipped:
20,040,354 bp masked across all 375 pack references that have sequence (0.348%), of which
104,184 bp is cloning vector in 45 contigs and 81,196 bp is cross-group sequence in 50.
The decision is recorded as D5, and D34 is why the denominator is now the whole pack:
twelve references had no masking result at all. **What is NOT done is the part D13
measured:** masking is detected at
`min_hit_len=200` / `min_identity=0.95` while reads map at 150 bp on shorter and less
identical stretches, so it cannot see read-scale recruitment - only 2 of 375 genomes carry
any host masking at all (2,737 bp), while the *T. gondii* reference absorbs ~5.4% of host
reads. Read-scale masking parameters are an untried lever with a real sensitivity cost; see
D13 before pulling it.

### B. Adopt UKHSA mSCAPE's control vocabulary and spike-in screen
mSCAPE (the most transparent public-health implementation; pipeline
`artic-network/scylla`, filters `ukhsa-collaboration/gpha-mscape-orangebox-claspar`)
carries an explicit metadata enum set:
  input_type          = community_standard | negative_control | positive_control |
                        specimen | validation_material
  control_type_details= NIBSC_11/242, NIBSC_20/170, bacillus_ms2phage,
                        resp_matrix_mc110, water_extraction_control, zepto_rp2.1,
                        zymo-mc_D6300
  spike_in            = ERCC-RNA_4456740, ms2-phage, phix, tobacco_mosaic_virus,
                        zymo_D6320/D6321, none
  spike_in_result     = pass | partial | fail
Our samplesheet has only role=sample|control, which cannot distinguish an extraction
blank from a no-template control from a mock community - and those fail differently.

**A `control_group` column was designed, measured and rejected for this slot (D21,
2026-10-13).** It would have scoped `control_rpm` per group so one scenario's
contaminated control could not suppress another's. Do not redesign it from scratch:
the measurement that killed it is that scoping `control_rpm` away from the whole run
REINSTATES D14 - re-scoring lod_v4 without its cross-scenario controls returns 30 rows,
15 of them as false `confirmed_in_silico` on E. coli, because `negative_host_only`'s
matrix blank deliberately relies on the run-global reach. Whatever `control_type`
becomes, it must leave a blank that covers the whole run expressible, and the two
mechanisms must be reconciled rather than stacked. See also D22: the `control_ratio`
multiplier that turns any control into a floor is itself underived.

### C. Negative-control-derived contaminant watchlist
Taxa present in >=50% of negative controls at >0.5 CPM go on a per-site watchlist,
**flagged not deleted** (decontam; Ibanez-Lligona et al. 2026, Front Microbiol,
doi:10.3389/fmicb.2026.1786985). decontam guidance is 5-6 negative controls. mSCAPE
ships an entire repo (`gpha-mscape-contamination`) just for aggregating negative
controls per site. We subtract against controls in the same run only.

### D. Machine-readable watchlist fetchers - only three sources are automatable
  ePATHogen  JSON/CSV API, refreshed nightly  -> the LEGAL source of Canadian Risk Groups
  EPPO       categorization.csv + v2 REST API, annual September revision
  eCFR API   US select agents
Everything else (Canada's notifiable list, CFIA animal and plant pest lists, WOAH
Chapter 1.3, WHO priority lists) is scrape-and-diff with no API.

**CORRECTION carried into the docs: HPTA Schedules 2-5 were REPEALED on 2026-03-26.**
The ePATHogen registry, not the Act's schedules, is now the legal source of Risk Group
assignments. Any text citing HPTA schedules as a watchlist source is out of date.

### E. HCID gate schema and mandatory human hand-off
mSCAPE's `hcid.json` uses `min_count`, `required_refs`, `min_coverage`,
`additional_refs`, `alt_taxon_id`, version-controlled with CI schema validation, and
their gate emits **reads plus per-reference `accession:reads:fraction_covered` rather
than a result**, with human hand-off mandatory. Our tiering is close in spirit; the
per-reference evidence emission and the explicit hand-off step are missing.

### F. Validation route for an UNTARGETED assay
  * WOAH Terrestrial Manual ch. 1.1.6 states the traditional validation pathway does
    not apply to untargeted assays and gives a QC-metric regime instead - the clearest
    official statement available for animal health.
  * EPPO has the most mature untargeted-method guidance of any body (plant health).
  * ENNGS/ESCV multi-site design: low-biomass-mimicking mock panel, multiple
    protocols, one central pipeline, LoD per organism, exclusivity testing against
    near neighbours AND organisms causing similar syndromes.
  * There is NO CLSI standard for metagenomic/untargeted sequencing (verified by
    enumeration of the MM series; MM24 2nd ed. 2021 covers targeted infectious-organism
    sequencing). No accredited mNGS EQA/proficiency scheme was found.
  * Precedent exists for an in-silico proficiency test of the dry lab alone
    (Front Microbiol 2020;11:575377) - directly applicable to our harness.

### F2. Divergent-detection layers, in the order the evidence supports
Run each as a RECALL-adding filter and push precision downstream into triage:
  1. Assemble, but **do not gate detection on assembly** - call at read level and
     assemble candidate reads for characterisation.
  2. Contig level: geNomad (sensitive) + Cenote-Taker 3, shown complementary. Do NOT
     take the union of five or six tools - measurably worse (mSystems 2024).
  3. Protein: DIAMOND --very-sensitive against RVDB-prot/NR, **plus an iterative
     MMseqs2 profile pass** (2-3 iterations): ROC1 0.41 -> 0.67 for remote homologs.
  4. RNA viruses: replace any Pfam-only RdRp screen with a **multi-database pHMM
     sweep at E <= 1e-5** (RdRpCATCH): F1 0.72 -> 0.996.
  5. Residual dark matter: ESMFold + Foldseek triage.
  6. Quality gate: CheckV 1.1.1.
Protein language models (LucaProt, Rider, LucaVirus) are competitive on RNA-virus
recall but **every published figure is the developers' own and no independent
head-to-head exists**; a 2026 preprint also shows pLM embedding space models viral
proteins worse than its general benchmarks imply. Not a standalone decision-maker.

### F3. ONT specifics
  * Require **>= 100 bp alignment length** for a taxonomic call (implemented).
  * Basecalling: sup@v5.2.0; never `fast`. There is no DNA sup v6.
  * Assemblers: myloasm or nanoMDBG over Flye for metagenomes.
  * Adopt **post-ligation pooling** if multiplexing low-biomass samples (barcode
    crosstalk).
  * Note: no 2025-2026 independent benchmark of ONT shotgun classifiers on real
    R10.4.1 data exists. sylph, KrakenUniq, Centrifuger and MetaPhlAn have
    essentially no independent ONT evidence.

### F4. Host depletion - reorder to protect-then-scrub
Move from single-reference Hostile to a pangenome or two-reference cascade. This is
the only published mitigation for losing integrated or host-homologous virus reads
(HPV, HBV, EBV, HERVs), and multi-reference filtering also blocks re-identification.
GRCh38 vs T2T-CHM13 for depletion is **contested** in the 2025 literature - Forbes
says T2T wins, the detaxizer benchmark says GRCh38. Unresolved, so do not assert it.

### G. Gaps where no published method exists
No general framework for vaccine-vs-wild-type calling from metagenomic data. No
published *Y. pestis* vs *Y. pseudotuberculosis* or *F. tularensis* subspecies
discrimination from metagenomes. No published reporting thresholds for select-agent
calls from metagenomes. No contamination-control guidance written for biothreat
metagenomics. Our marker tier is addressing real gaps, which also means it is
unvalidated and must not be presented as established practice.

## Watchlist provenance (do before operational use)
The seed watchlist is derived from published notifiable/high-consequence lists but is
NOT authoritative. Reconcile against, and record the version of:
- Human: PHAC ePATHogen risk-group database; national notifiable disease list
- Animal: WOAH listed diseases
- Plant: EPPO A1/A2

Each entry's `source` column should name the list and its version once reconciled.
