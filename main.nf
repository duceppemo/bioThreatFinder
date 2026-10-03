#!/usr/bin/env nextflow
/*
 * bioThreatFinder
 *
 * Two-stage metagenomic detection of high-consequence pathogens.
 *
 *   Stage 1  broad k-mer / ANI screen  -> candidate taxa (permissive)
 *   Stage 2  targeted competitive alignment -> evidence (strict)
 *
 * The split exists because stage-1 classifiers are built for community
 * profiling, not for single-organism calls. A handful of k-mers shared with an
 * environmental relative is enough to put a high-consequence taxon in a Kraken2
 * report. That is a lead, not a detection. Stage 2 is what decides.
 */

nextflow.enable.dsl = 2

include { FASTP; CHOPPER; SEQKIT_STATS      } from './modules/local/preprocess'
include { MINIMAP2_HOST_DEPLETE             } from './modules/local/preprocess'
include { RESOLVE_WATCHLIST                 } from './modules/local/detect'
include { KRAKEN2; BRACKEN; SYLPH_PROFILE   } from './modules/local/detect'
include { KRAKENUNIQ; METABULI              } from './modules/local/detect'
include { PARSE_SCREEN                      } from './modules/local/detect'
include { MINIMAP2_COMPETITIVE              } from './modules/local/detect'
include { COVERAGE_STATS                    } from './modules/local/detect'
include { MARKER_ALIGN_NT; MARKER_SEARCH_AA } from './modules/local/detect'
include { MARKER_SEARCH_HMM; MARKER_EVIDENCE } from './modules/local/detect'
include { SCORE_EVIDENCE; MAKE_REPORT       } from './modules/local/detect'

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

def helpMessage() {
    log.info """
    bioThreatFinder v${workflow.manifest.version}

    Usage:
      nextflow run . -profile conda,workstation \\
          --input samplesheet.csv \\
          --pack human \\
          --outdir results

    Required:
      --input            Samplesheet CSV: sample,platform,fastq_1,fastq_2,role
                         platform = illumina | ont
                         role     = sample | control   (controls drive contamination subtraction)
      --pack             Domain pack: human | animal | plant   [${params.pack}]

    Key options:
      --outdir           Output directory                      [${params.outdir}]
      --db_root          Shared read-only DBs (taxdump, sylph, nt)  [${params.db_root}]
      --btf_db           Project-owned DBs on fast NVMe             [${params.btf_db}]

    Stage-1 classifier arms (any-classifier OR; consensus recorded, not required):
      --skip_kraken2     Disable the PRIMARY arm. Kraken2 with
                         --report-minimizer-data supplies distinct-minimizer counts
                         and a duplicity ratio, which is what separates reads spread
                         across a genome from reads piled on a conserved locus.
      --skip_krakenuniq false
                         OPT-IN second opinion (default off - note that OMITTING this
                         flag leaves the arm OFF; `--krakenuniq` is not a parameter and
                         Nextflow accepts unknown params silently). Kraken2's minimizer
                         columns are literally KrakenUniq's HyperLogLog code at the
                         same default precision, counting minimizers rather than
                         k-mers, and the only published comparison finds comparable
                         F1 once the threshold is shifted. KrakenUniq has had no
                         commit since Dec 2023. Prefer the prebuilt index over a
                         build if you enable this.
      --skip_metabuli    Disable the DNA+amino-acid arm. Independent benchmarks
                         rank its precision poorly and several upstream bugs are
                         open, so treat its calls as leads to confirm.
      --skip_sylph       Disable ANI/containment.
      --min_distinct_kmers / --max_kmer_dup   stage-1 k-mer quality flags
      --skip_host_depletion
      --skip_markers     Disable marker gating. Taxa whose calls depend on
                         plasmid/toxin loci will then be promoted on alignment
                         evidence alone, which is NOT defensible for several.
      Stage-2 threshold OVERRIDES. packs/<pack>/pack.yaml is the source of truth;
      these only override it, and null means "use the pack value". An unknown key is
      rejected rather than ignored.
      --min_unique_reads    on-target reads      [${params.min_unique_reads ?: 'pack.yaml'}]
      --min_breadth         reference breadth    [${params.min_breadth ?: 'pack.yaml'}]
      --max_depth_cv        depth CV ceiling     [${params.max_depth_cv ?: 'pack.yaml'}]
      --min_read_identity   identity floor       [${params.min_read_identity ?: 'pack.yaml (per platform AND group)'}]

    Profiles: conda, mamba, docker, singularity, workstation, test
    """.stripIndent()
}

/* Pack YAML carries ${db_root} placeholders so a pack is portable between
 * machines without editing it. Resolve them here, once. */
def resolvePlaceholders(Object v, Map ctx) {
    if (v instanceof String) {
        def out = v
        ctx.each { k, val -> out = out.replace('${' + k + '}', val.toString()) }
        return out
    }
    if (v instanceof Map)  return v.collectEntries { k, val -> [k, resolvePlaceholders(val, ctx)] }
    if (v instanceof List) return v.collect { resolvePlaceholders(it, ctx) }
    return v
}

def loadPack(String packName) {
    def packDir = file("${params.packs_dir}/${packName}")
    if (!packDir.exists()) {
        exit 1, "Pack '${packName}' not found at ${packDir}. Available: " +
                file(params.packs_dir).list().sort().join(', ')
    }
    def yml = packDir.resolve('pack.yaml')
    if (!yml.exists()) exit 1, "Pack '${packName}' has no pack.yaml (${yml})"
    def parsed = new org.yaml.snakeyaml.Yaml().load(yml.text)
    def resolved = resolvePlaceholders(
        parsed,
        [db_root: params.db_root, btf_db: params.btf_db, projectDir: projectDir])
    resolved.dir = packDir
    return resolved
}

/* Pads an optional per-sample channel out to cover every sample.
 *
 * `base` carries (meta, placeholder) for all samples; `produced` carries
 * (meta, file) only for samples where the arm actually ran. A remainder join keeps
 * the missing ones and substitutes the placeholder, so a skipped or failed
 * classifier degrades that one arm rather than dropping the sample from the screen
 * entirely - which would look like a clean negative. */
def fillMissing(base, produced) {
    return base
        .join(produced, remainder: true)
        .map { meta, placeholder, actual -> tuple(meta, actual ?: placeholder) }
}

/* Guards against the worst failure mode this pipeline has.
 *
 * Classifier output is joined to packs/<pack>/watchlist.resolved.tsv, which is keyed
 * by NCBI taxid. Metabuli publishes GTDB-based indexes as well, and GTDB taxids are
 * NOT NCBI taxids. Point the pipeline at a GTDB index and every watchlist lookup
 * misses - producing not an error but an empty result, indistinguishable from a
 * clean negative across the entire screen.
 *
 * A silent false negative on every watchlist organism is worse than a crash, so this
 * is refused rather than warned about. The build scripts record the namespace in
 * BTF_BUILD_MANIFEST.txt precisely so it can be checked here. */
def assertNcbiTaxonomy(dbDir, label) {
    def manifest = file("${dbDir}/BTF_BUILD_MANIFEST.txt")
    if (!manifest.exists()) {
        log.warn "${label}: no BTF_BUILD_MANIFEST.txt under ${dbDir}, so the taxonomy " +
                 "namespace cannot be verified. If this index is GTDB-based, its taxids " +
                 "will not join the NCBI-keyed watchlist and this arm will silently " +
                 "contribute nothing. Build it with scripts/build_metabuli_db.sh to get " +
                 "a manifest, or pass --skip_metabuli to be explicit."
        return
    }
    def line = manifest.readLines().find { it.startsWith('taxonomy_namespace=') }
    if (!line) {
        log.warn "${label}: BTF_BUILD_MANIFEST.txt under ${dbDir} has no " +
                 "taxonomy_namespace entry; namespace could not be verified."
        return
    }
    def ns = line.split('=', 2)[1].trim()
    def nsUpper = ns.toUpperCase()

    // "NCBI" ALONE IS NOT ENOUGH OF AN ANSWER.
    //
    // build_metabuli_db.sh used to derive this line from --mode, so a GTDB prebuilt
    // unpacked under --mode refseq was stamped NCBI and sailed through here. The
    // builder now checks the written taxonomy instead, and says when it COULD NOT
    // verify ("NCBI (UNVERIFIED ...)") or when the artifact and the request disagree
    // ("NCBI (CONTRADICTION ...)"). Both still start with NCBI, so accepting any
    // NCBI-prefixed string would re-open the hole one layer up.
    if (nsUpper.contains('UNVERIFIED') || nsUpper.contains('CONTRADICTION')) {
        exit 1, """
        ERROR: the ${label} index at ${dbDir} records an UNCONFIRMED taxonomy namespace:

          ${ns}

        A GTDB-namespaced index joined to an NCBI-keyed watchlist yields ZERO matches -
        not an error, just an empty result indistinguishable from a clean negative for
        every organism on the watchlist. An unverified namespace is not evidence that
        the index is usable, so it is refused rather than assumed.

        Resolve by one of:
          - rebuild so the taxonomy can be checked: scripts/build_metabuli_db.sh --mode refseq
          - point --metabuli_db at an index whose manifest verifies as NCBI
          - drop the arm deliberately:              --skip_metabuli
        """.stripIndent()
    }
    if (!nsUpper.startsWith('NCBI')) {
        exit 1, """
        ERROR: the ${label} index at ${dbDir} is built in the '${ns}' taxonomy namespace.

        This pipeline joins classifier output to the watchlist by NCBI taxid. A
        non-NCBI index yields ZERO watchlist matches - not an error, just an empty
        result, which is indistinguishable from a clean negative for every organism
        on the watchlist. That is refused rather than warned about.

        Resolve by one of:
          - build an NCBI-namespaced index:  scripts/build_metabuli_db.sh --mode refseq
          - point --metabuli_db at an NCBI-namespaced index
          - drop the arm deliberately:       --skip_metabuli
        """.stripIndent()
    }
}

def parseSamplesheet(path) {
    return Channel.fromPath(path, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def required = ['sample', 'platform', 'fastq_1']
            required.each { c ->
                if (!row.containsKey(c)) exit 1, "Samplesheet is missing column '${c}'"
            }
            def id = row.sample?.trim()
            if (!id) exit 1, "Samplesheet has a row with an empty 'sample' value"

            def platform = (row.platform ?: '').trim().toLowerCase()
            if (!(platform in ['illumina', 'ont'])) {
                exit 1, "Sample '${id}': platform must be 'illumina' or 'ont', got '${row.platform}'"
            }

            def role = ((row.role ?: 'sample').trim().toLowerCase()) ?: 'sample'
            if (!(role in ['sample', 'control'])) {
                exit 1, "Sample '${id}': role must be 'sample' or 'control', got '${row.role}'"
            }

            def r1 = file(row.fastq_1.trim(), checkIfExists: true)
            def r2 = (row.fastq_2?.trim()) ? file(row.fastq_2.trim(), checkIfExists: true) : null

            if (platform == 'ont' && r2) {
                exit 1, "Sample '${id}': platform 'ont' cannot have fastq_2"
            }

            def meta = [id: id, platform: platform, role: role, single_end: (r2 == null)]
            return tuple(meta, r2 ? [r1, r2] : [r1])
        }
        // A HEADER-ONLY SAMPLESHEET MUST NOT LOOK LIKE A SUCCESSFUL RUN.
        //
        // tests/data/samplesheet_test.csv was one line - the header - so `-profile
        // test` yielded no samples: every per-sample channel was empty, SCORE_EVIDENCE
        // ran on .ifEmpty([]) inputs, MAKE_REPORT emitted an empty report and the run
        // exited 0. The profile documented as existing "so the channel wiring can be
        // exercised in minutes" exercised none of it, and a green result was
        // indistinguishable from a working pipeline.
        .ifEmpty { exit 1, """
        ERROR: samplesheet ${params.input} contains no data rows.

        Nothing would be analysed and the run would still exit 0, which is
        indistinguishable from a clean negative for every organism on the watchlist.
        Add at least one sample row: sample,platform,fastq_1,fastq_2,role
        """.stripIndent() }
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    if (params.help) { helpMessage(); exit 0 }
    if (!params.input) { helpMessage(); exit 1, "ERROR: --input is required" }

    def pack = loadPack(params.pack)

    def arms = []
    if (!params.skip_krakenuniq) arms << 'krakenuniq'
    if (!params.skip_metabuli)   arms << 'metabuli'
    if (!params.skip_kraken2)    arms << 'kraken2/bracken'
    if (!params.skip_sylph)      arms << 'sylph'
    if (arms.isEmpty()) {
        exit 1, "ERROR: every stage-1 classifier is disabled. Nothing would be screened."
    }
    if (params.skip_kraken2) {
        log.warn "Kraken2 is disabled. It is the PRIMARY nucleotide arm and the source " +
                 "of distinct-minimizer evidence (--report-minimizer-data), so without " +
                 "it reads piled on a few conserved k-mers cannot be told apart from " +
                 "reads spread across a genome at stage 1."
    }

    log.info """
    bioThreatFinder v${workflow.manifest.version}
      pack            : ${pack.name}
      stage-1 arms    : ${arms.join(', ')}  (any-classifier OR; consensus recorded, not required)
      marker gating   : ${params.skip_markers ? 'DISABLED' : 'enabled'}
      outdir          : ${params.outdir}
      report          : ${params.outdir}/report/biothreatfinder_report.html
      findings TSV    : ${params.outdir}/score_evidence/findings.tsv

    Findings are in-silico screening evidence only. Every one requires orthogonal
    laboratory confirmation and interpretation by a qualified laboratorian before
    any clinical or public-health action.
    """.stripIndent()

    ch_reads = parseSamplesheet(params.input)

    // -- Watchlist resolution --------------------------------------------
    // Runs once per pipeline invocation, not per sample. Fails the whole run if
    // any watchlist name is unresolvable, before any compute is spent.
    RESOLVE_WATCHLIST(
        file("${pack.dir}/${pack.reporting.watchlist}", checkIfExists: true),
        file(params.taxdump, checkIfExists: true)
    )
    ch_watchlist   = RESOLVE_WATCHLIST.out.resolved
    ch_descendants = RESOLVE_WATCHLIST.out.descendants

    // -- Stage 0: QC ------------------------------------------------------
    ch_reads.branch { meta, _reads ->
        illumina: meta.platform == 'illumina'
        ont:      meta.platform == 'ont'
    }.set { ch_by_platform }

    FASTP(ch_by_platform.illumina)
    CHOPPER(ch_by_platform.ont)

    ch_trimmed = FASTP.out.reads.mix(CHOPPER.out.reads)

    // -- Stage 0b: host depletion -----------------------------------------
    // Host reads are the majority of most clinical/veterinary specimens and a
    // major source of spurious classification via repeats.
    // A run with --skip_host_depletion must not inherit the previous run's QC file.
    // storeDir writes outside the work dir and is not cleaned between runs, so without a
    // sentinel, run B (skipped) would present run A's healthy depletion figures as its own.
    if (params.skip_host_depletion) {
        ch_clean = ch_trimmed
        ch_host_stats = Channel.empty()
        Channel.of("sample\tplatform\talignment_records\thost_matching_records\t" +
                   "host_fraction_of_records\tsurviving_records\n" +
                   "ALL_SAMPLES\tNA\tNA\tNA\tSKIPPED\tNA\n")
            .collectFile(name: 'host_depletion.tsv', storeDir: "${params.outdir}/qc")
            .subscribe { f -> log.warn "[QC] host depletion SKIPPED; eukaryotic references " +
                               "recruit host reads in proportion to what survives it " +
                               "(docs/decisions.md D14) -> ${f}" }
    } else {
        ch_host_idx = ch_trimmed.map { meta, reads ->
            def idx = meta.platform == 'illumina'
                    ? pack.host.minimap2_index_sr
                    : pack.host.minimap2_index_lr
            tuple(meta, reads, file(idx, checkIfExists: true))
        }
        MINIMAP2_HOST_DEPLETE(ch_host_idx)
        ch_clean      = MINIMAP2_HOST_DEPLETE.out.reads
        ch_host_stats = MINIMAP2_HOST_DEPLETE.out.stats
        // Run-level QC artifact. Previously ch_host_stats was produced and consumed by
        // nothing; D14 showed the number is load-bearing, so it is now written out.
        MINIMAP2_HOST_DEPLETE.out.efficiency
            .collectFile(name: 'host_depletion.tsv',
                         storeDir: "${params.outdir}/qc",
                         keepHeader: true, skip: 1, sort: true)
            .subscribe { f -> log.info "[QC] host depletion -> ${f}" }
    }

    SEQKIT_STATS(ch_clean)

    // -- Stage 1: broad screen (multi-classifier) --------------------------
    // Independent arms, because they fail in different directions. Kraken2 is
    // PRIMARY (decisions.md D2); KrakenUniq is opt-in and OFF by default:
    //   kraken2     minimizer counts + duplicity -> specificity at stage 1, primary
    //   metabuli    DNA + amino acid             -> divergent / novel agents
    //   bracken     abundance re-estimation      -> abundance, continuity
    //   sylph       ANI / containment            -> orthogonal, low FP
    //   krakenuniq  HyperLogLog cardinality      -> OPT-IN second opinion
    // Any arm raises a candidate; consensus is recorded, not required. Requiring
    // consensus would discard the sensitivity that running four arms buys.
    // AN ABSENT OPTIONAL PATH IS AN EMPTY LIST, NOT A PLACEHOLDER FILE.
    //
    // This used to be `def no_file = file("assets/NO_FILE")`, and every padded arm got
    // that same file. Nextflow stages inputs by basename, so a process receiving two or
    // more of them died with "input file name collision" - which is what PARSE_SCREEN did
    // the moment two stage-1 arms were skipped. stageAs papered over it by giving each a
    // private directory; passing `[]` removes the cause, because nothing is staged at all
    // and there is no basename to collide. It is also the nf-core idiom, and the script
    // side becomes `x ? ... : ''` instead of comparing a filename to a magic string.

    // pack.yaml supplies the database paths; a CLI --*_db overrides it.
    def ku_db_path = params.krakenuniq_db ?: pack.screen.krakenuniq_db
    def mb_db_path = params.metabuli_db   ?: pack.screen.metabuli_db
    def k2_db_path = params.kraken2_db    ?: pack.screen.kraken2_db
    def sy_db_path = params.sylph_db      ?: pack.screen.sylph_db

    // Every arm is padded out to all samples by fillMissing(), so a skipped or
    // failed classifier never silently drops a sample from the join below.
    ch_arm_base = ch_clean.map { meta, _r -> tuple(meta, []) }

    if (params.skip_krakenuniq) {
        ch_ku = ch_clean.map { meta, _r -> tuple(meta, []) }
    } else {
        KRAKENUNIQ(ch_clean, file(ku_db_path, checkIfExists: true))
        ch_ku = fillMissing(ch_arm_base, KRAKENUNIQ.out.report)
    }

    if (params.skip_metabuli) {
        ch_mb = ch_clean.map { meta, _r -> tuple(meta, []) }
    } else {
        def mb_db = file(mb_db_path, checkIfExists: true)
        // Refuses a GTDB-namespaced index outright - see assertNcbiTaxonomy().
        assertNcbiTaxonomy(mb_db, 'Metabuli')
        METABULI(ch_clean, mb_db)
        ch_mb = fillMissing(ch_arm_base, METABULI.out.report)
    }

    if (params.skip_kraken2) {
        ch_k2      = ch_clean.map { meta, _r -> tuple(meta, []) }
        ch_bracken = ch_clean.map { meta, _r -> tuple(meta, []) }
    } else {
        KRAKEN2(ch_clean, file(k2_db_path, checkIfExists: true))
        ch_k2 = fillMissing(ch_arm_base, KRAKEN2.out.report)
        // Bracken is built on a fixed read length, so it applies to short reads
        // only; ONT read lengths violate its model.
        BRACKEN(KRAKEN2.out.report.filter { meta, _r -> meta.platform == 'illumina' },
                file(k2_db_path))
        ch_bracken = fillMissing(ch_arm_base, BRACKEN.out.report)
    }

    if (params.skip_sylph) {
        ch_sylph = ch_clean.map { meta, _r -> tuple(meta, []) }
    } else {
        SYLPH_PROFILE(ch_clean, file(sy_db_path, checkIfExists: true))
        ch_sylph = fillMissing(ch_arm_base, SYLPH_PROFILE.out.profile)
    }

    ch_screen = ch_ku.join(ch_k2).join(ch_bracken).join(ch_mb).join(ch_sylph)

    PARSE_SCREEN(ch_screen, ch_watchlist.first(), ch_descendants.first())

    // -- Stage 2: targeted confirmation -----------------------------------
    // Fan out one alignment job per (sample, candidate taxon). Candidates are
    // few by construction, so this stays cheap even on large runs.
    // NO SECOND READ-COUNT GATE HERE. THIS ONE SILENTLY DISCARDED EVERY SYLPH HIT.
    //
    // There used to be a
    //     .filter { (row.assigned_reads as Integer) >= params.min_screen_reads }
    // which looks like a harmless echo of the stage-1 threshold and is neither
    // harmless nor an echo. `assigned_reads` is a READ COUNT, and only the read-based
    // arms (kraken2, krakenuniq, metabuli) produce one; sylph reports ANI containment
    // and abundance, so a sylph-raised candidate carries assigned_reads=0 BY
    // CONSTRUCTION. With min_screen_reads=10, `0 >= 10` is false and the candidate was
    // dropped before stage 2 ever saw it.
    //
    // Measured on the first end-to-end run: Y. pestis at sylph abundance 6.46 and
    // E. coli at 22.17 were both raised correctly by PARSE_SCREEN, both filtered out
    // here, MINIMAP2_COMPETITIVE received nothing, and the run exited 0 with a
    // header-only findings.tsv - on a positive control whose scenario file says "any
    // outcome other than confirmed is a configuration fault". A clean negative for the
    // entire watchlist, reported as success.
    //
    // bin/parse_screen.py is the single place candidacy is decided: it applies
    // --min-reads to each read-based classifier, requires non-empty support, and
    // documents "Permissive OR: any classifier raises the candidate. Stage 2 decides."
    // Re-filtering on a read count here overrode that and contradicted it.
    ch_candidates = PARSE_SCREEN.out.candidates
        .splitCsv(header: true, sep: '\t', elem: 1)
        .map { meta, row -> tuple(meta, row) }

    ch_confirm_in = ch_candidates
        .combine(ch_clean, by: 0)
        .map { meta, row, reads ->
            // Competitive reference bundle: target + near neighbours, built by
            // scripts/build_pack_refs.py.
            //
            // Prefer the prebuilt, preset-matched minimap2 index. k and w are
            // baked in at index time, so an -x sr index performs badly on ONT
            // reads and vice versa - hence one index per preset rather than one
            // shared. Fall back to the FASTA (minimap2 indexes on the fly) so a
            // refs dir built with --skip-index still works.
            def preset = meta.platform == 'illumina' ? 'sr' : 'map-ont'
            def mmi = file("${pack.confirm.refs_dir}/competitive/${row.taxid}.${preset}.mmi")
            def ref = mmi.exists()
                    ? mmi
                    : file("${pack.confirm.refs_dir}/competitive/${row.taxid}.fna.gz")
            tuple(meta + [taxid: row.taxid, taxon: row.name, group: row.group,
                          assigned_reads: row.assigned_reads], reads, ref, row)
        }
        .branch { _meta, _reads, ref, _row ->
            present: ref.exists()
            missing: true
        }

    // A watchlist taxon with no reference genome cannot be confirmed. That is a
    // gap in the database, not a negative result, and it has to reach the report.
    ch_missing_ref = ch_confirm_in.missing.map { meta, _reads, ref, row ->
        tuple(meta, row + [confirm_status: 'no_reference', reference_path: ref.toString()])
    }

    MINIMAP2_COMPETITIVE(ch_confirm_in.present.map { meta, reads, ref, _row -> tuple(meta, reads, ref) })
    // Soft-mask intervals, when the pack has been masked. Absent is tolerated and
    // reported rather than assumed clean - an unmasked reference set is a known
    // false-positive source, not a neutral default.
    def masked_dir = pack.confirm.masked_dir ? file(pack.confirm.masked_dir) : []
    if (!masked_dir || !masked_dir.exists()) {
        log.warn "No masked reference intervals at ${pack.confirm.masked_dir}. " +
                 "References are UNMASKED, and vector, host, low-complexity or " +
                 "cross-contaminated sequence inside a reference will collect reads " +
                 "and inflate breadth. Run scripts/mask_pack_refs.py --pack ${pack.name}."
        masked_dir = []
    }
    COVERAGE_STATS(MINIMAP2_COMPETITIVE.out.bam, masked_dir)

    // -- Stage 2b: marker tier --------------------------------------------
    // Run per sample, not per candidate: the marker DB is small, so one pass
    // serves every taxon's rule more cheaply than a fan-out.
    if (params.skip_markers) {
        ch_marker_states = Channel.empty()
    } else {
        ch_markers_dir = file(pack.confirm.markers_dir ?: "${pack.confirm.refs_dir}/../markers")
        ch_marker_cat  = file("${pack.dir}/markers.tsv", checkIfExists: true)

        MARKER_ALIGN_NT(ch_clean, ch_markers_dir)
        MARKER_SEARCH_AA(ch_clean, ch_markers_dir)

        if (params.skip_marker_hmm) {
            ch_hmm = ch_clean.map { meta, _r -> tuple(meta, []) }
        } else {
            MARKER_SEARCH_HMM(MARKER_SEARCH_AA.out.hit_reads, ch_markers_dir)
            ch_hmm = MARKER_SEARCH_HMM.out.hits
        }

        ch_marker_in = MARKER_ALIGN_NT.out.bam
            .join(MARKER_SEARCH_AA.out.hits, remainder: true)
            .join(ch_hmm, remainder: true)
            .map { meta, bam, bai, aa, hmm ->
                tuple(meta, bam ?: [], bai ?: [], aa ?: [], hmm ?: [])
            }

        MARKER_EVIDENCE(ch_marker_in, ch_marker_cat, ch_markers_dir)
        ch_marker_states = MARKER_EVIDENCE.out.states
    }

    // -- Evidence scoring, control subtraction, tiering -------------------
    SCORE_EVIDENCE(
        COVERAGE_STATS.out.stats.map { meta, tsv -> tsv }.collect().ifEmpty([]),
        PARSE_SCREEN.out.candidates.map { meta, tsv -> tsv }.collect().ifEmpty([]),
        ch_missing_ref.map { meta, row -> "${meta.id}\t${meta.taxid}\t${meta.taxon}\tno_reference" }
            .collectFile(name: 'missing_refs.tsv', newLine: true).ifEmpty([]),
        ch_watchlist.first(),
        file("${pack.dir}/pack.yaml"),
        ch_marker_states.map { meta, tsv -> tsv }.collect().ifEmpty([]),
        params.skip_markers ? [] : file("${pack.dir}/marker_rules.yaml", checkIfExists: true)
    )

    MAKE_REPORT(
        SCORE_EVIDENCE.out.findings,
        SEQKIT_STATS.out.stats.map { meta, s -> s }.collect().ifEmpty([]),
        file("${pack.dir}/pack.yaml")
    )
}

// NOTE: the top-level `workflow.onComplete` handler this file used to carry is not
// permitted by the Nextflow 24+ strict script syntax, which allows only
// declarations at script level. The run notice now prints from inside the entry
// workflow instead (the run-summary block in the entry workflow), and the caveats live where
// they actually matter - on every finding in the HTML report and the findings TSV,
// not in terminal scrollback that nobody keeps.
