/*
 * QC and host depletion. Platform-specific by design: short and long reads need
 * different trimmers and different minimap2 presets, and conflating them is a
 * reliable way to produce garbage.
 */

process FASTP {
    tag   "${meta.id}"
    label 'process_medium'
    conda "bioconda::fastp=1.3.7"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.trim*.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.fastp.json")    , emit: json
    tuple val(meta), path("${meta.id}.fastp.html")    , emit: html
    path "versions.yml"                               , emit: versions

    script:
    // Adapter and low-complexity trimming matter more here than in a typical
    // assembly workflow: adapter dimers and poly-G tails are a documented source
    // of spurious taxonomic assignment.
    def io = meta.single_end
        ? "-i ${reads[0]} -o ${meta.id}.trim.fastq.gz"
        : "-i ${reads[0]} -I ${reads[1]} -o ${meta.id}.trim_1.fastq.gz -O ${meta.id}.trim_2.fastq.gz --detect_adapter_for_pe"
    """
    fastp \\
        ${io} \\
        --length_required ${params.sr_min_len} \\
        --qualified_quality_phred ${params.sr_min_qual} \\
        --low_complexity_filter \\
        --trim_poly_g \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        2> ${meta.id}.fastp.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
    END_VERSIONS
    """
}

process CHOPPER {
    tag   "${meta.id}"
    label 'process_medium'
    conda "bioconda::chopper=0.14.1"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.trim.fastq.gz"), emit: reads
    path "versions.yml"                              , emit: versions

    script:
    """
    # Nanopore quality filtering. --headcrop removes the adapter-proximal bases
    # where basecalling accuracy is lowest.
    gzip -cd ${reads} \\
      | chopper --quality ${params.lr_min_qual} --minlength ${params.lr_min_len} --headcrop 10 --threads ${task.cpus} \\
      | gzip -c > ${meta.id}.trim.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: \$(chopper --version 2>&1 | sed 's/chopper //')
    END_VERSIONS
    """
}

process MINIMAP2_HOST_DEPLETE {
    tag   "${meta.id}"
    label 'process_high'
    // python is here for bin/host_depletion_stats.py, which turns flagstat into the
    // host-depletion figure D14 shows is load-bearing. Same pin as MAKE_REPORT.
    conda "bioconda::minimap2=2.31 bioconda::samtools=1.22 conda-forge::python=3.11"

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path("${meta.id}.dehost*.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.host_stats.txt")  , emit: stats
    path "${meta.id}.host_depletion.tsv"                , emit: efficiency
    path "versions.yml"                                 , emit: versions

    script:
    def preset = meta.platform == 'illumina' ? 'sr' : 'map-ont'
    // Paired end: -f 12 keeps only pairs where BOTH mates are unmapped. Keeping
    // half-mapped pairs would drag host sequence into classification through the
    // unmapped mate.
    def extract = meta.single_end
        ? "samtools fastq -f 4 -0 ${meta.id}.dehost.fastq.gz -"
        : "samtools fastq -f 12 -F 256 -1 ${meta.id}.dehost_1.fastq.gz -2 ${meta.id}.dehost_2.fastq.gz -0 /dev/null -s /dev/null -"
    """
    minimap2 -ax ${preset} -t ${task.cpus} --secondary=no ${index} ${reads} 2> ${meta.id}.minimap2.log \\
      | tee >(samtools flagstat - > ${meta.id}.host_stats.txt) \\
      | ${extract}

    echo "# host index: ${index}" >> ${meta.id}.host_stats.txt

    # Host-depletion efficiency as a first-class artifact. D14: eukaryotic references
    # recruit host reads in proportion to what survives this step - the T. gondii
    # reference takes 7.74% of them - so this fraction gates a whole class of false
    # positive. flagstat already computed it, host_stats was published, and ch_host_stats
    # was a dead channel that nothing read.
    host_depletion_stats.py \\
        --flagstat ${meta.id}.host_stats.txt \\
        --sample '${meta.id}' \\
        --platform '${meta.platform}' \\
        --out ${meta.id}.host_depletion.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}

process SEQKIT_STATS {
    tag   "${meta.id}"
    label 'process_low'
    conda "bioconda::seqkit=2.14.0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.seqkit_stats.tsv"), emit: stats
    path "versions.yml"                                 , emit: versions

    script:
    """
    seqkit stats -a -T -j ${task.cpus} ${reads} > ${meta.id}.seqkit_stats.tsv
    # Prepend sample identity so downstream aggregation does not have to parse
    # filenames back into sample names.
    sed -i "1s/\$/\\tsample\\trole\\tplatform/" ${meta.id}.seqkit_stats.tsv
    awk 'NR>1{print \$0"\\t${meta.id}\\t${meta.role}\\t${meta.platform}"}' ${meta.id}.seqkit_stats.tsv > tmp \\
      && head -1 ${meta.id}.seqkit_stats.tsv > final && cat tmp >> final && mv final ${meta.id}.seqkit_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """
}
