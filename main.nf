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
      --db_root          Shared database root                   [${params.db_root}]
      --skip_host_depletion / --skip_sylph
      --min_unique_reads Stage-2 minimum on-target reads        [${params.min_unique_reads}]
      --min_breadth      Stage-2 minimum reference breadth      [${params.min_breadth}]

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
    def resolved = resolvePlaceholders(parsed, [db_root: params.db_root, projectDir: projectDir])
    resolved.dir = packDir
    return resolved
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
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    if (params.help) { helpMessage(); exit 0 }
    if (!params.input) { helpMessage(); exit 1, "ERROR: --input is required" }

    def pack = loadPack(params.pack)
    log.info "bioThreatFinder v${workflow.manifest.version}  pack=${pack.name}  outdir=${params.outdir}"

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
    if (params.skip_host_depletion) {
        ch_clean = ch_trimmed
        ch_host_stats = Channel.empty()
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
    }

    SEQKIT_STATS(ch_clean)

    // -- Stage 1: broad screen --------------------------------------------
    KRAKEN2(ch_clean, file(pack.screen.kraken2_db, checkIfExists: true))

    // Bracken re-estimates abundance from Kraken2's k-mer assignments. It is
    // built on a fixed read length, so it applies to short reads only; ONT read
    // lengths violate its model and it is skipped there.
    ch_for_bracken = KRAKEN2.out.report.filter { meta, _r -> meta.platform == 'illumina' }
    BRACKEN(ch_for_bracken, file(pack.screen.kraken2_db))

    if (params.skip_sylph) {
        ch_sylph = Channel.empty()
    } else {
        // Orthogonal evidence: sylph is containment/ANI based, so it disagrees
        // with Kraken2 in exactly the cases that matter - a taxon called from
        // scattered conserved k-mers gets no sylph support.
        SYLPH_PROFILE(ch_clean, file(pack.screen.sylph_db, checkIfExists: true))
        ch_sylph = SYLPH_PROFILE.out.profile
    }

    // Join the screens per sample. Sylph and Bracken are both optional, so pad
    // with a placeholder rather than dropping the sample from the join.
    ch_screen = KRAKEN2.out.report
        .join(KRAKEN2.out.classified_out)
        .join(BRACKEN.out.report, remainder: true)
        .join(ch_sylph, remainder: true)
        .map { meta, report, kout, bracken, sylph ->
            tuple(meta, report, kout,
                  bracken ?: file("${projectDir}/assets/NO_FILE"),
                  sylph   ?: file("${projectDir}/assets/NO_FILE"))
        }

    PARSE_SCREEN(ch_screen, ch_watchlist.first(), ch_descendants.first())

    // -- Stage 2: targeted confirmation -----------------------------------
    // Fan out one alignment job per (sample, candidate taxon). Candidates are
    // few by construction, so this stays cheap even on large runs.
    ch_candidates = PARSE_SCREEN.out.candidates
        .splitCsv(header: true, sep: '\t', elem: 1)
        .map { meta, row -> tuple(meta, row) }
        .filter { _meta, row -> (row.assigned_reads as Integer) >= params.min_screen_reads }

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
            tuple(meta + [taxid: row.taxid, taxon: row.name, group: row.group], reads, ref, row)
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
    COVERAGE_STATS(MINIMAP2_COMPETITIVE.out.bam)

    // -- Stage 2b: marker tier --------------------------------------------
    // Run per sample, not per candidate: the marker DB is small, so one pass
    // serves every taxon's rule more cheaply than a fan-out.
    def no_file = file("${projectDir}/assets/NO_FILE")
    if (params.skip_markers) {
        ch_marker_states = Channel.empty()
    } else {
        ch_markers_dir = file(pack.confirm.markers_dir ?: "${pack.confirm.refs_dir}/../markers")
        ch_marker_cat  = file("${pack.dir}/markers.tsv", checkIfExists: true)

        MARKER_ALIGN_NT(ch_clean, ch_markers_dir)
        MARKER_SEARCH_AA(ch_clean, ch_markers_dir)

        if (params.skip_marker_hmm) {
            ch_hmm = ch_clean.map { meta, _r -> tuple(meta, no_file) }
        } else {
            MARKER_SEARCH_HMM(MARKER_SEARCH_AA.out.hit_reads, ch_markers_dir)
            ch_hmm = MARKER_SEARCH_HMM.out.hits
        }

        ch_marker_in = MARKER_ALIGN_NT.out.bam
            .join(MARKER_SEARCH_AA.out.hits, remainder: true)
            .join(ch_hmm, remainder: true)
            .map { meta, bam, bai, aa, hmm ->
                tuple(meta, bam ?: no_file, bai ?: no_file,
                      aa ?: no_file, hmm ?: no_file)
            }

        MARKER_EVIDENCE(ch_marker_in, ch_marker_cat, ch_markers_dir)
        ch_marker_states = MARKER_EVIDENCE.out.states
    }

    // -- Evidence scoring, control subtraction, tiering -------------------
    SCORE_EVIDENCE(
        COVERAGE_STATS.out.stats.map { meta, tsv -> tsv }.collect().ifEmpty([]),
        PARSE_SCREEN.out.candidates.map { meta, tsv -> tsv }.collect().ifEmpty([]),
        ch_missing_ref.map { meta, row -> "${meta.id}\t${meta.taxid}\t${meta.taxon}\tno_reference" }
            .collectFile(name: 'missing_refs.tsv', newLine: true).ifEmpty(file("${projectDir}/assets/NO_FILE")),
        ch_watchlist.first(),
        file("${pack.dir}/pack.yaml"),
        ch_marker_states.map { meta, tsv -> tsv }.collect().ifEmpty([]),
        params.skip_markers ? no_file : file("${pack.dir}/marker_rules.yaml", checkIfExists: true)
    )

    MAKE_REPORT(
        SCORE_EVIDENCE.out.findings,
        SEQKIT_STATS.out.stats.map { meta, s -> s }.collect().ifEmpty([]),
        file("${pack.dir}/pack.yaml")
    )
}

workflow.onComplete {
    log.info """
    ${workflow.success ? 'Completed' : 'FAILED'}  duration=${workflow.duration}
    Report: ${params.outdir}/report/biothreatfinder_report.html
    Findings TSV: ${params.outdir}/score_evidence/findings.tsv

    Reminder: findings are in-silico screening evidence only and require
    orthogonal laboratory confirmation before any clinical or public-health action.
    """.stripIndent()
}
