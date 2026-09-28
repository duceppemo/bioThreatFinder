# Roadmap

## Blocking a real run
1. **Upgrade Nextflow** — 21.10.4 installed, >=23.10 required.
2. **Create the tool env**
   `mamba create -n biothreatfinder -c bioconda -c conda-forge kraken2 bracken minimap2 samtools seqkit fastp chopper sylph ncbi-datasets-cli`
   `datasets`, `kraken2-build`, `bracken-build` and `seqkit` are all absent from this box.
   The only `minimap2` is the one bundled with ont-dorado-server.
3. **Build databases** — `scripts/fetch_host_genome.sh`, `build_pack_refs.py`, `build_kraken2_db.sh`.
   Kraken2 PlusPF needs ~700 GiB transient against 1074 GiB free. Tight but passes.
4. **Validate on known-positive data** before trusting any output. Mock communities
   (ZymoBIOMICS) plus spiked dilution series to establish an actual limit of detection.
   Until this is done the thresholds in `packs/*/pack.yaml` are reasoned guesses, not
   measured operating points.

## Next feature, in priority order

### DONE - Marker-gene / plasmid module
Implemented in v0.2.0: 163 markers, 58 rule sets, three-valued rule engine,
three evidence tiers (dna/protein/hmm). Every high_consequence taxon has a rule.
See README for the design and validation. Remaining marker work:
  * 31 notifiable taxa still have no rule
  * curate_required markers: B. mallei absent regions, vaccinia-specific loci,
    F. novicida discriminating loci, C. auris clades
  * absence-region logic as a first-class marker class (B. mallei needs it)

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

### 3. Accuracy upgrades to stages 1 and 2
  * Stage 1: KrakenUniq (unique k-mer counts via HyperLogLog) or Metabuli
    (joint DNA + amino acid) instead of Kraken2 alone - cuts species-level false
    positives at the source
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

## Watchlist provenance (do before operational use)
The seed watchlist is derived from published notifiable/high-consequence lists but is
NOT authoritative. Reconcile against, and record the version of:
- Human: PHAC ePATHogen risk-group database; national notifiable disease list
- Animal: WOAH listed diseases
- Plant: EPPO A1/A2

Each entry's `source` column should name the list and its version once reconciled.
