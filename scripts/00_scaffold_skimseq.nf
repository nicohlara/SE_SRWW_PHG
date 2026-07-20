nextflow.enable.dsl=2

// baseline parameters
params.fastq_dir = "/90daydata/guedira_seq_map/skim_seq"
params.output_dir = "/90daydata/guedira_seq_map/nico/ragtag_assemblies"

//files
params.ref = "/90daydata/guedira_seq_map/RefCS_2.1/iwgsc_refseqv2.1_assembly.fa"

workflow {
    Channel
        .fromPath("${params.fastq_dir}/*_.fastq.gz")
        .map { fastq ->
            tuple(fastq.baseName.replaceFirst(/_interleaved$/, ''), fastq)
        }
        .set { fastq_ch }

    create_contigs(fastq_ch) |
        scaffold
}

process create_contigs {
    input:
    tuple val(sample_id), path(fastq)
	
    output:
    tuple val(sample_id), path("${sample_id}_contigs.fasta")

    script:
    """
    ##assemble contigs from raw reads
    hifiasm \
        -o ${sample_id} \
        -f 0 \
        -l 0 \
        ${fastq}

    awk '/^S/{print ">"\$2; print \$3}' ${sample_id}.asm.bp.hap1.p_ctg.gfa > ${sample_id}_contigs.fasta
    """
}

process scaffold {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)

    output:
    path "${sample_id}_scaffold"

    script:
    """
    ragtag.py scaffold \
        -q 60 \
        -f 100000 \
        -i 0.6 \
        --remove-small \
        -o ${sample_id}_scaffold \
        ${params.ref} \
        ${fastq}
    """
}
