/*
 * Stage 1 (screen), stage 2 (confirm), scoring and reporting.
 */

process RESOLVE_WATCHLIST {
    label 'process_low'
    conda "conda-forge::python=3.11"
    publishDir "${params.outdir}/watchlist", mode: params.publish_dir_mode

    input:
    path watchlist
    path taxdump

    output:
    path "watchlist.resolved.tsv"   , emit: resolved
    path "watchlist.descendants.tsv", emit: descendants
    path "resolve.log"              , emit: log

    script:
    // No --allow-unresolved here on purpose. An unresolvable watchlist name means
    // an organism is silently not being screened for, which is the worst possible
    // failure mode for this tool. Fail the run before spending any compute.
    """
    resolve_watchlist.py \\
        --watchlist ${watchlist} \\
        --taxdump ${taxdump} \\
        --require-fresh-taxdump \\
        --max-taxdump-age-days ${params.max_taxdump_age_days} \\
        --out watchlist.resolved.tsv \\
        --descendants-out watchlist.descendants.tsv \\
        2> resolve.log
    cat resolve.log >&2
    """
}

process KRAKEN2 {
    tag   "${meta.id}"
    label 'process_hugemem'
    conda "bioconda::kraken2=2.17.2"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("${meta.id}.kraken2.report.txt"), emit: report
    tuple val(meta), path("${meta.id}.kraken2.output.txt"), emit: classified_out
    path "versions.yml"                                   , emit: versions

    script:
    def paired = meta.single_end ? '' : '--paired'
    """
    kraken2 \\
        --db ${db} \\
        --threads ${task.cpus} \\
        ${paired} \\
        --gzip-compressed \\
        --report ${meta.id}.kraken2.report.txt \\
        --output ${meta.id}.kraken2.output.txt \\
        --report-minimizer-data \\
        ${reads} \\
        2> ${meta.id}.kraken2.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(kraken2 --version | head -1 | sed 's/Kraken version //')
    END_VERSIONS
    """
}

process BRACKEN {
    tag   "${meta.id}"
    label 'process_low'
    conda "bioconda::bracken=3.1p1"

    // Bracken fails outright when a report has no reads at the requested level.
    // That is a legitimate outcome for a clean negative control, so it must not
    // kill the run.
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(report)
    path db

    output:
    tuple val(meta), path("${meta.id}.bracken.report.txt"), emit: report
    path "versions.yml"                                   , emit: versions

    script:
    """
    # --report-minimizer-data adds columns Bracken cannot parse; strip them.
    awk 'BEGIN{FS=OFS="\\t"} {if (NF>6) {\$4=""; \$5=""; gsub(/\\t\\t/,"\\t")} print}' ${report} > k2_for_bracken.txt || cp ${report} k2_for_bracken.txt

    bracken \\
        -d ${db} \\
        -i k2_for_bracken.txt \\
        -o ${meta.id}.bracken.tsv \\
        -w ${meta.id}.bracken.report.txt \\
        -r ${params.bracken_readlen ?: 150} \\
        -l S \\
        -t 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: \$(bracken -v 2>&1 | sed 's/Bracken v//')
    END_VERSIONS
    """
}

process SYLPH_PROFILE {
    tag   "${meta.id}"
    label 'process_high'
    conda "bioconda::sylph=1.0.0"

    errorStrategy 'ignore'

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("${meta.id}.sylph.tsv"), emit: profile
    path "versions.yml"                          , emit: versions

    script:
    def input_args = meta.single_end ? "-r ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    # -u gives full containment-based abundance rather than the default
    # coverage-adjusted estimate; we want the raw ANI/containment evidence.
    sylph profile ${db} ${input_args} -t ${task.cpus} -o ${meta.id}.sylph.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph --version | sed 's/sylph //')
    END_VERSIONS
    """
}

/*
 * ---------------------------------------------------------------------------
 * Stage 1, OPT-IN arm: KrakenUniq (--skip_krakenuniq defaults true)
 *
 * Reports distinct K-MER cardinality per taxon (HyperLogLog) plus a duplication
 * ratio and estimated coverage. Kraken2's --report-minimizer-data columns are
 * literally this same HyperLogLog implementation at the same default precision,
 * counting distinct MINIMIZERS - a subsample - which is why Kraken2 is primary and
 * this is a second opinion rather than a dependency: the only published comparison
 * finds comparable F1 once the threshold is shifted by the 2-3x that subsampling
 * implies, and KrakenUniq has had no commit since December 2023.
 *
 * If you do enable it, download a prebuilt index rather than building one: a build
 * costs ~1400 GiB transient and needs Jellyfish 1.x. See docs/decisions.md D2.
 * ---------------------------------------------------------------------------
 */
process KRAKENUNIQ {
    tag   "${meta.id}"
    label 'process_krakenuniq'
    conda "bioconda::krakenuniq=1.0.4"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("${meta.id}.krakenuniq.report.txt"), emit: report
    tuple val(meta), path("${meta.id}.krakenuniq.output.txt.gz"), emit: classified_out
    path "versions.yml"                                       , emit: versions

    script:
    def paired = meta.single_end ? '' : '--paired'
    // Loading the index once is far faster per sample than memory-mapping it
    // repeatedly, which matters across a benchmark batch. The flag name differs
    // between KrakenUniq versions, so it is a configurable string rather than a
    // guess baked in here - see params.krakenuniq_preload_args.
    def preload = params.krakenuniq_preload_args ?: ''
    def hll_override = params.krakenuniq_hll_precision ?: ''
    """
    # HyperLogLog precision is a CLASSIFY-time parameter, not a build parameter:
    # the sketches are built while reading each sample. build_krakenuniq_db.sh
    # records the intended value beside the index, so query and build cannot
    # drift. An explicit --krakenuniq_hll_precision overrides it.
    HLL_ARGS=""
    if [ -n "${hll_override}" ]; then
        HLL_ARGS="--hll-precision ${hll_override}"
    elif [ -f "${db}/btf_classify_defaults.env" ]; then
        # shellcheck disable=SC1091
        . "${db}/btf_classify_defaults.env"
        # KRAKENUNIQ_HLL_PRECISION is the key build_krakenuniq_db.sh writes into
        # btf_classify_defaults.env. This tested `hll_precision_for_classify`, which
        # appears only inside BTF_BUILD_MANIFEST.txt - a file this module never reads -
        # so HLL_ARGS was ALWAYS empty and krakenuniq fell back to its own default.
        # Building with --hll-precision 14 and classifying at 12 is exactly the
        # build/query drift this mechanism exists to prevent, and it was silent.
        # Latent only because the build default (12) happens to match krakenuniq's.
        if [ -n "\${KRAKENUNIQ_HLL_PRECISION:-}" ]; then
            HLL_ARGS="--hll-precision \${KRAKENUNIQ_HLL_PRECISION}"
        fi
    fi
    echo "krakenuniq HLL args: \${HLL_ARGS:-<tool default>}" >&2

    krakenuniq \
        --db ${db} \
        --threads ${task.cpus} \
        ${paired} \
        ${preload} \
        \${HLL_ARGS} \
        --report-file ${meta.id}.krakenuniq.report.txt \
        --output ${meta.id}.krakenuniq.output.txt \
        ${reads} \
        2> ${meta.id}.krakenuniq.log

    # The per-read output is large and only needed for read extraction; keep it
    # compressed rather than discarding it, so an analyst can pull the reads
    # behind a finding later.
    gzip -f ${meta.id}.krakenuniq.output.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakenuniq: \$(krakenuniq --version 2>&1 | head -1 | sed 's/^.*version //')
    END_VERSIONS
    """
}

/*
 * ---------------------------------------------------------------------------
 * Stage 1, divergence arm: Metabuli
 *
 * Joint DNA + amino-acid k-mer classification. Nucleotide k-mers miss a divergent
 * genome entirely - a novel orthopoxvirus or arenavirus at 75% nucleotide identity
 * to anything in the database is essentially invisible to Kraken2/KrakenUniq, but
 * conserved in protein space. Until this arm existed, only the marker tier could
 * see through divergence, and only for curated markers.
 * ---------------------------------------------------------------------------
 */
process METABULI {
    tag   "${meta.id}"
    label 'process_metabuli'
    conda "bioconda::metabuli=1.2.0"

    // A Metabuli failure must not sink the run: Kraken2 is the primary arm (D2) and
    // the screen parser records which classifiers actually produced output.
    //
    // CAVEAT on that last clause: screen_summary.tsv is the record, and it reaches
    // neither SCORE_EVIDENCE nor MAKE_REPORT - it is only published. For a taxon that
    // WAS raised, n_classifiers/classifiers do travel into findings.tsv, so a missing
    // arm is visible there. But if this arm were the ONLY one that would have raised a
    // taxon, its ignored failure loses that taxon with no trace in the report at all.
    // The summary is published beside the results as the compensating record.
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("${meta.id}.metabuli.report.tsv"), emit: report
    path "versions.yml"                                    , emit: versions

    script:
    // VERIFY: --seq-mode values and the positional argument order against the
    // installed Metabuli version. 1 = single-end, 2 = paired-end, 3 = long read.
    def seq_mode = meta.single_end ? (meta.platform == 'ont' ? '3' : '1') : '2'
    def input_args = meta.single_end ? "${reads}" : "${reads[0]} ${reads[1]}"
    """
    mkdir -p mb_out
    metabuli classify \
        --seq-mode ${seq_mode} \
        --threads ${task.cpus} \
        ${input_args} \
        ${db} \
        mb_out \
        ${meta.id} \
        2> ${meta.id}.metabuli.log

    # VERIFY: output filename convention (<jobid>_report.tsv).
    if [ -s mb_out/${meta.id}_report.tsv ]; then
        cp mb_out/${meta.id}_report.tsv ${meta.id}.metabuli.report.tsv
    else
        # Produce an empty report rather than failing the task, so the join
        # downstream still has a file and the parser records metabuli as absent.
        : > ${meta.id}.metabuli.report.tsv
        echo "WARNING: metabuli produced no report" >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabuli: \$(metabuli 2>&1 | grep -i version | head -1 | sed 's/^.*ersion[: ]*//')
    END_VERSIONS
    """
}

process PARSE_SCREEN {
    tag   "${meta.id}"
    label 'process_low'
    conda "conda-forge::python=3.11"
    publishDir "${params.outdir}/screen", mode: params.publish_dir_mode

    // EACH ARM STAGES INTO ITS OWN DIRECTORY.
    //
    // fillMissing() pads every absent arm with the SAME assets/NO_FILE placeholder, and
    // Nextflow stages inputs by basename - so with two or more arms skipped the task got
    // several inputs called NO_FILE and died with "input file name collision". Running
    // with --skip_kraken2 --skip_metabuli (krakenuniq off by default) makes four of
    // these NO_FILE at once. It stayed latent under the default configuration, where
    // only krakenuniq is absent on Illumina, which is exactly why no review caught it:
    // the padding mechanism that exists to stop a sample being silently dropped broke
    // the moment more than one arm was missing.
    //
    // stageAs also protects against two REAL reports sharing a basename, which the
    // per-arm filenames do not currently guarantee.
    input:
    tuple val(meta),
          path(ku_report, stageAs: 'arm_krakenuniq/*'),
          path(k2_report, stageAs: 'arm_kraken2/*'),
          path(bracken,   stageAs: 'arm_bracken/*'),
          path(metabuli,  stageAs: 'arm_metabuli/*'),
          path(sylph,     stageAs: 'arm_sylph/*')
    path watchlist
    path descendants

    output:
    tuple val(meta), path("${meta.id}.candidates.tsv"), emit: candidates
    tuple val(meta), path("${meta.id}.screen_summary.tsv"), optional: true, emit: summary
    path "${meta.id}.screen.log"                      , emit: log

    script:
    def ku_arg      = ku_report.name != 'NO_FILE' ? "--krakenuniq-report ${ku_report}" : ''
    def k2_arg      = k2_report.name != 'NO_FILE' ? "--kraken2-report ${k2_report}"    : ''
    def bracken_arg = bracken.name   != 'NO_FILE' ? "--bracken-report ${bracken}"      : ''
    def mb_arg      = metabuli.name  != 'NO_FILE' ? "--metabuli-report ${metabuli}"     : ''
    def sylph_arg   = sylph.name     != 'NO_FILE' ? "--sylph-profile ${sylph}"          : ''
    """
    parse_screen.py \\
        ${ku_arg} \\
        ${k2_arg} \\
        ${bracken_arg} \\
        ${mb_arg} \\
        ${sylph_arg} \\
        --min-distinct-kmers ${params.min_distinct_kmers} \\
        --max-kmer-dup ${params.max_kmer_dup} \\
        --watchlist ${watchlist} \\
        --descendants ${descendants} \\
        --sample ${meta.id} \\
        --role ${meta.role} \\
        --platform ${meta.platform} \\
        --min-reads ${params.min_screen_reads} \\
        --out ${meta.id}.candidates.tsv \\
        --out-summary ${meta.id}.screen_summary.tsv \\
        2> ${meta.id}.screen.log
    cat ${meta.id}.screen.log >&2
    """
}

process MINIMAP2_COMPETITIVE {
    tag   "${meta.id}:${meta.taxid}"
    label 'process_medium'
    conda "bioconda::minimap2=2.31 bioconda::samtools=1.22"

    input:
    tuple val(meta), path(reads), path(ref)

    output:
    tuple val(meta), path("${meta.id}.${meta.taxid}.bam"),
                     path("${meta.id}.${meta.taxid}.bam.bai"),
                     path("${meta.id}.${meta.taxid}.nreads"), emit: bam
    path "versions.yml"                                     , emit: versions

    script:
    def preset = meta.platform == 'illumina' ? 'sr' : 'map-ont'
    """
    # Competitive mapping against target + near neighbours in one index. Reads
    # that fit a neighbour better go to the neighbour, and reads that fit both
    # equally get MAPQ 0 and are discarded downstream. --secondary=no keeps one
    # primary alignment per read so the assignment is unambiguous.
    minimap2 -ax ${preset} -t ${task.cpus} --secondary=no ${ref} ${reads} 2> mm2.log \\
      | samtools sort -@ ${task.cpus} -o ${meta.id}.${meta.taxid}.bam -
    samtools index -@ ${task.cpus} ${meta.id}.${meta.taxid}.bam

    # Denominator for reads-per-million, so a finding can be compared across
    # samples of different depth and against the negative control.
    total=0
    for f in ${reads}; do
        n=\$( { gzip -cd "\$f" 2>/dev/null || cat "\$f"; } | awk 'END{print int(NR/4)}' )
        total=\$(( total + n ))
    done
    echo "\$total" > ${meta.id}.${meta.taxid}.nreads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}

process COVERAGE_STATS {
    tag   "${meta.id}:${meta.taxid}"
    label 'process_low'
    conda "bioconda::samtools=1.22 conda-forge::python=3.11"
    publishDir "${params.outdir}/confirm", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai), path(nreads)
    path masked_dir

    output:
    tuple val(meta), path("${meta.id}.${meta.taxid}.coverage.tsv"), emit: stats

    script:
    def mask_arg = masked_dir.name != 'NO_FILE' ? "--mask-bed ${masked_dir}" : ''
    """
    coverage_stats.py \\
        --bam ${bam} \\
        --sample ${meta.id} \\
        --taxid ${meta.taxid} \\
        --taxon '${meta.taxon}' \\
        --group '${meta.group}' \\
        --role ${meta.role} \\
        --platform ${meta.platform} \\
        --mapq-min ${params.mapq_min} \\
        --total-reads \$(cat ${nreads}) \\
        --assigned-reads ${meta.assigned_reads ?: 0} \\
        ${mask_arg} \\
        --out ${meta.id}.${meta.taxid}.coverage.tsv
    """
}

process SCORE_EVIDENCE {
    label 'process_low'
    conda "conda-forge::python=3.11"
    publishDir "${params.outdir}/score_evidence", mode: params.publish_dir_mode

    // Same collision hazard as PARSE_SCREEN: missing_refs and marker_rules are both
    // NO_FILE when there are no database gaps AND marker gating is off, which would
    // stage two files called NO_FILE into one directory.
    input:
    path coverage
    path candidates
    path missing_refs, stageAs: 'missing_refs/*'
    path watchlist
    path pack_yaml
    path marker_states
    path marker_rules, stageAs: 'marker_rules/*' 

    output:
    path "findings.tsv", emit: findings
    path "score.log"   , emit: log

    script:
    def cov_arg  = coverage   ? "--coverage ${coverage.join(' ')}"     : ''
    def cand_arg = candidates ? "--candidates ${candidates.join(' ')}" : ''
    def miss_arg = missing_refs.name != 'NO_FILE' ? "--missing-refs ${missing_refs}" : ''
    def mk_arg   = marker_states ? "--marker-evidence ${marker_states.join(' ')}" : ''
    def mr_arg   = marker_rules.name != 'NO_FILE' ? "--marker-rules ${marker_rules}" : ''
    // pack.yaml is the source of truth for thresholds; these params exist so a
    // profile can relax a gate without editing the pack. They were DECLARED in
    // nextflow.config and passed to nothing, so `-profile test` asserted relaxed
    // gates that never applied and its empty result looked like a clean negative.
    // Null means "use the pack value". score_evidence.py rejects an unknown key.
    def thr = [
        'min_unique_reads': params.min_unique_reads,
        'min_breadth': params.min_breadth,
        'max_depth_cv': params.max_depth_cv,
        'min_read_identity': params.min_read_identity,
        'mapq_min': params.mapq_min,
    ].findAll { _k, v -> v != null }
     .collect { k, v -> "--threshold ${k}=${v}" }
     .join(' ')
    """
    score_evidence.py \\
        ${cov_arg} \\
        ${cand_arg} \\
        ${miss_arg} \\
        ${mk_arg} \\
        ${mr_arg} \\
        --watchlist ${watchlist} \\
        --pack-yaml ${pack_yaml} \\
        --control-ratio ${params.control_ratio} \\
        ${thr} \\
        --out findings.tsv \\
        2> score.log
    cat score.log >&2
    """
}

process MAKE_REPORT {
    label 'process_low'
    conda "conda-forge::python=3.11"
    publishDir "${params.outdir}/report", mode: params.publish_dir_mode

    input:
    path findings
    path read_stats
    path pack_yaml

    output:
    path "biothreatfinder_report.html", emit: html
    path "biothreatfinder_summary.tsv", emit: summary

    script:
    def stats_arg = read_stats ? "--read-stats ${read_stats.join(' ')}" : ''
    """
    make_report.py \\
        --findings ${findings} \\
        ${stats_arg} \\
        --pack-yaml ${pack_yaml} \\
        --pack ${params.pack} \\
        --lab '${params.lab_name}' \\
        --run-id '${params.run_id ?: workflow.runName}' \\
        --pipeline-version '${workflow.manifest.version}' \\
        --command '${workflow.commandLine.replace("'", "")}' \\
        --html biothreatfinder_report.html \\
        --summary biothreatfinder_summary.tsv
    """
}

/*
 * ---------------------------------------------------------------------------
 * Marker tier (stage 2b)
 *
 * Runs once per sample against the combined marker database rather than once
 * per candidate taxon: the marker DB is small, so one pass is cheaper than a
 * fan-out and gives every taxon's rule the evidence it needs.
 *
 * Three tiers, complementary rather than a consensus vote - see
 * bin/marker_evidence.py for why.
 * ---------------------------------------------------------------------------
 */

process MARKER_ALIGN_NT {
    tag   "${meta.id}"
    label 'process_medium'
    conda "bioconda::minimap2=2.31 bioconda::samtools=1.22"

    input:
    tuple val(meta), path(reads)
    path markers_dir

    output:
    tuple val(meta), path("${meta.id}.markers.nt.bam"), path("${meta.id}.markers.nt.bam.bai"), emit: bam

    script:
    def preset = meta.platform == 'illumina' ? 'sr' : 'map-ont'
    """
    # Markers are short (single genes / regions), so --secondary=no plus the
    # MAPQ filter downstream is what keeps a paralogue from claiming a hit.
    minimap2 -ax ${preset} -t ${task.cpus} --secondary=no \\
        ${markers_dir}/markers.nt.fna ${reads} 2> mm2.log \\
      | samtools sort -@ ${task.cpus} -o ${meta.id}.markers.nt.bam -
    samtools index -@ ${task.cpus} ${meta.id}.markers.nt.bam
    """
}

process MARKER_SEARCH_AA {
    tag   "${meta.id}"
    label 'process_high'
    conda "bioconda::mmseqs2=18.8cc5c bioconda::seqkit=2.14.0"

    input:
    tuple val(meta), path(reads)
    path markers_dir

    output:
    tuple val(meta), path("${meta.id}.markers.aa.m8"), emit: hits
    tuple val(meta), path("${meta.id}.marker_hit_reads.fasta.gz"), emit: hit_reads

    script:
    """
    cat ${reads} > all.fastq.gz
    seqkit fq2fa all.fastq.gz -o reads.fasta

    # Translated search. This is the tier that sees through synonymous
    # divergence, so an allele that nucleotide alignment misses is still found.
    # -s 7.5 is near-maximum sensitivity; on an HPC node the cost is acceptable
    # and a missed divergent toxin gene is not.
    # THE HMM TIER MUST SEE MORE THAN THE PROTEIN TIER ALREADY FOUND.
    #
    # hit_reads used to be built from this search at -e 1e-5, i.e. only reads the
    # protein tier had ALREADY called a hit on. MARKER_SEARCH_HMM therefore ran on a
    # subset of what the protein tier found and could only ever confirm a subset of it -
    # so a read divergent enough to NEED a profile search never reached hmmsearch, and
    # the marker came back ABSENT with no indication the HMM tier had never seen it.
    # The README describes this tier as "the reason a novel BoNT or Stx subtype, or a
    # divergent filovirus gene, is still caught", which as wired it structurally could
    # not do.
    #
    # One search at the LOOSE cutoff, then two views of it: rows at or below
    # marker_aa_evalue are the protein tier's hits (unchanged semantics), while the
    # candidate read set for the HMM tier comes from every row. Same compute as before -
    # the cost is in the search, not the cutoff.
    mmseqs easy-search reads.fasta ${markers_dir}/markers.aa.faa \\
        all_hits.m8 tmp_mmseqs \\
        --threads ${task.cpus} \\
        -s 7.5 \\
        -e ${params.marker_hmm_prefilter_evalue} \\
        --format-output 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits' \\
        2> mmseqs.log

    # The protein tier's own hits: the strict cutoff, applied to the shared search.
    awk -F'\\t' -v e=${params.marker_aa_evalue} \\
        'NF >= 11 && (\$11 + 0) <= e' all_hits.m8 > ${meta.id}.markers.aa.m8 || true

    # Candidate reads for the HMM tier: every row of the loose search, so a read too
    # divergent to clear the protein cutoff still gets a profile search. Still a small
    # subset of the library - six-frame translating every read would cost far more than
    # it returns.
    cut -f1 all_hits.m8 | sort -u > hit_ids.txt || true
    if [ -s hit_ids.txt ]; then
        seqkit grep -f hit_ids.txt reads.fasta | gzip -c > ${meta.id}.marker_hit_reads.fasta.gz
    else
        : | gzip -c > ${meta.id}.marker_hit_reads.fasta.gz
    fi
    rm -rf tmp_mmseqs reads.fasta all.fastq.gz all_hits.m8
    """
}

process MARKER_SEARCH_HMM {
    tag   "${meta.id}"
    label 'process_medium'
    conda "bioconda::hmmer=3.4 bioconda::seqkit=2.14.0"

    // An empty hit set is a legitimate outcome, not a failure.
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(hit_reads)
    path markers_dir

    output:
    tuple val(meta), path("${meta.id}.markers.hmm.domtbl"), emit: hits

    script:
    """
    touch ${meta.id}.markers.hmm.domtbl
    if [ -s ${hit_reads} ] && [ "\$(gzip -cd ${hit_reads} | head -c1 | wc -c)" != "0" ]; then
        # Six-frame translation: read orientation and frame are unknown, and
        # missing a frame would silently lose the hit the HMM tier exists for.
        seqkit translate -f 6 --allow-unknown-codon ${hit_reads} -o translated.faa 2>> seqkit.log
        if [ -s translated.faa ]; then
            hmmsearch --cpu ${task.cpus} \\
                --domtblout ${meta.id}.markers.hmm.domtbl \\
                -E 1e-3 --domE 1e-3 \\
                ${markers_dir}/markers.hmm translated.faa > hmmsearch.log 2>&1 || true
        fi
    fi
    """
}

process MARKER_EVIDENCE {
    tag   "${meta.id}"
    label 'process_low'
    conda "bioconda::samtools=1.22 conda-forge::python=3.11"
    publishDir "${params.outdir}/markers", mode: params.publish_dir_mode

    // aa_hits and hmm_hits are both NO_FILE when the protein and HMM tiers are skipped
    // or their (errorStrategy 'ignore') processes failed - two inputs, one basename.
    input:
    tuple val(meta),
          path(nt_bam,   stageAs: 'tier_nt/*'),
          path(nt_bai,   stageAs: 'tier_nt/*'),
          path(aa_hits,  stageAs: 'tier_aa/*'),
          path(hmm_hits, stageAs: 'tier_hmm/*')
    path markers_catalogue
    path markers_dir

    output:
    tuple val(meta), path("${meta.id}.marker_states.tsv"), emit: states
    path "${meta.id}.marker.log"                         , emit: log

    script:
    def nt_arg  = nt_bam.name   != 'NO_FILE' ? "--nt-bam ${nt_bam}"     : ''
    def aa_arg  = aa_hits.name  != 'NO_FILE' ? "--aa-hits ${aa_hits} --aa-fasta ${markers_dir}/markers.aa.faa" : ''
    def hmm_arg = hmm_hits.name != 'NO_FILE' ? "--hmm-hits ${hmm_hits}" : ''
    """
    marker_evidence.py \\
        --markers ${markers_catalogue} \\
        --manifest ${markers_dir}/markers.manifest.tsv \\
        --sample ${meta.id} \\
        ${nt_arg} ${aa_arg} ${hmm_arg} \\
        --mapq-min ${params.mapq_min} \\
        --platform ${meta.platform} \\
        --hmm-min-bits ${params.hmm_min_bits} \\
        --min-low-depth-reads ${params.marker_min_low_depth_reads} \\
        --low-depth-breadth-fraction ${params.marker_low_depth_breadth_fraction} \\
        --out ${meta.id}.marker_states.tsv \\
        2> ${meta.id}.marker.log
    cat ${meta.id}.marker.log >&2
    """
}
