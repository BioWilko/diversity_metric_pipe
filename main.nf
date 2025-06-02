process align_to_ref {

    container "community.wave.seqera.io/library/bwa-mem2_minimap2_samtools:649dd658311f2f07"

    cpus 8
    memory '16 GB'
    time '1h'

    input:
    tuple val(sample), path(reads1), path(reads2), path(ref)

    output:
    tuple val(sample), path("${sample}.bam")

    script:
    preset_str = reads2 ? "-x sr" : "-x map-ont"

    index_ref_str = reads2 ? "bwa-mem2 index -p ${ref} ${ref}" : ""

    alignment_str = reads2 ? "bwa-mem2 mem -t ${task.cpus} ${ref} ${reads1} ${reads2}" : "minimap2 ${preset_str} -a -t ${task.cpus} ${ref} ${reads1}"

    """
    ${index_ref_str}
    ${alignment_str} | samtools view -bS -F 4 - | samtools sort -o ${sample}.bam
    """
}

process maptide {
    container "community.wave.seqera.io/library/pip_maptide_setuptools:8df043169de12124"

    cpus 1
    memory "16 GB"
    time "30m"

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.depths.tsv")

    script:
    """
    maptide ${bam} > ${sample}.depths.tsv
    """
}

process nucleotide_diversity {

    container "python:3.10.17-slim-bookworm"

    cpus 1
    memory "8GB"
    time "1h"

    input:
    tuple val(sample), path(depth_tsv)

    output:
    tuple val(sample), path("${sample}.diversity_metrics.csv")

    script:
    """
    nucleotide_diversity.py ${sample} ${depth_tsv} > ${sample}.diversity_metrics.csv
    """
}

workflow {
    Channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            [row.sample, row.reads1, row.reads2, row.ref]
        }
        .map { sample, reads1, reads2, ref ->
            [sample, reads1, reads2 ? reads2 : [], ref]
        }
        .set { samplesheet_ch }

    align_to_ref(samplesheet_ch)
        | maptide
        | nucleotide_diversity

    nucleotide_diversity.out
        .map { _sample, diversity_csv ->
            diversity_csv
        }
        .collectFile(name: "diversity_metrics.csv", storeDir: params.outdir, keepHeader: true, skip: 1)
}
