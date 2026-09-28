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
        --out watchlist.resolved.tsv \\
        --descendants-out watchlist.descendants.tsv \\
        2> resolve.log
    cat resolve.log >&2
    """
}

process KRAKEN2 {
    tag   "${meta.id}"
    label 'process_hugemem'
    conda "bioconda::kraken2=2.1.3"

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
    conda "bioconda::bracken=2.9"

    input:
    tuple val(meta), path(report)
    path db

    output:
    tuple val(meta), path("${meta.id}.bracken.report.txt"), emit: report
    path "versions.yml"                                   , emit: versions

    // Bracken fails outright when a report has no reads at the requested level.
    // That is a legitimate outcome for a clean negative control, so it must not
    // kill the run.
    errorStrategy 'ignore'

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
    conda "bioconda::sylph=0.6.1"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("${meta.id}.sylph.tsv"), emit: profile
    path "versions.yml"                          , emit: versions

    errorStrategy 'ignore'

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

process PARSE_SCREEN {
    tag   "${meta.id}"
    label 'process_low'
    conda "conda-forge::python=3.11"
    publishDir "${params.outdir}/screen", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(k2_report), path(k2_out), path(bracken), path(sylph)
    path watchlist
    path descendants

    output:
    tuple val(meta), path("${meta.id}.candidates.tsv"), emit: candidates
    path "${meta.id}.screen.log"                      , emit: log

    script:
    def bracken_arg = bracken.name != 'NO_FILE' ? "--bracken-report ${bracken}" : ''
    def sylph_arg   = sylph.name   != 'NO_FILE' ? "--sylph-profile ${sylph}"   : ''
    """
    parse_screen.py \\
        --kraken2-report ${k2_report} \\
        ${bracken_arg} \\
        ${sylph_arg} \\
        --watchlist ${watchlist} \\
        --descendants ${descendants} \\
        --sample ${meta.id} \\
        --role ${meta.role} \\
        --platform ${meta.platform} \\
        --min-reads ${params.min_screen_reads} \\
        --out ${meta.id}.candidates.tsv \\
        2> ${meta.id}.screen.log
    cat ${meta.id}.screen.log >&2
    """
}

process MINIMAP2_COMPETITIVE {
    tag   "${meta.id}:${meta.taxid}"
    label 'process_medium'
    conda "bioconda::minimap2=2.28 bioconda::samtools=1.20"

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
    conda "bioconda::samtools=1.20 conda-forge::python=3.11"
    publishDir "${params.outdir}/confirm", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai), path(nreads)

    output:
    tuple val(meta), path("${meta.id}.${meta.taxid}.coverage.tsv"), emit: stats

    script:
    """
    coverage_stats.py \\
        --bam ${bam} \\
        --sample ${meta.id} \\
        --taxid ${meta.taxid} \\
        --taxon '${meta.taxon}' \\
        --group '${meta.group}' \\
        --role ${meta.role} \\
        --mapq-min ${params.mapq_min} \\
        --total-reads \$(cat ${nreads}) \\
        --out ${meta.id}.${meta.taxid}.coverage.tsv
    """
}

process SCORE_EVIDENCE {
    label 'process_low'
    conda "conda-forge::python=3.11"
    publishDir "${params.outdir}/score_evidence", mode: params.publish_dir_mode

    input:
    path coverage
    path candidates
    path missing_refs
    path watchlist
    path pack_yaml

    output:
    path "findings.tsv", emit: findings
    path "score.log"   , emit: log

    script:
    def cov_arg  = coverage   ? "--coverage ${coverage.join(' ')}"     : ''
    def cand_arg = candidates ? "--candidates ${candidates.join(' ')}" : ''
    def miss_arg = missing_refs.name != 'NO_FILE' ? "--missing-refs ${missing_refs}" : ''
    """
    score_evidence.py \\
        ${cov_arg} \\
        ${cand_arg} \\
        ${miss_arg} \\
        --watchlist ${watchlist} \\
        --pack-yaml ${pack_yaml} \\
        --control-ratio ${params.control_ratio} \\
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
