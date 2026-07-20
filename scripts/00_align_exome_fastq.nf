nextflow.enable.dsl=2

// baseline parameters
params.fastq_dir = "/90daydata/guedira_seq_map/mwillman/exome_capture_fastq_data/merged_by_line"
//params.fastq_dir = "/90daydata/guedira_seq_map/nico/exome_test"
params.output_dir = "/90daydata/guedira_seq_map/nico/exome_capture_vcf"

//files
params.ref = "/90daydata/guedira_seq_map/RefCS_2.1/iwgsc_refseqv2.1_assembly.fa"

workflow {
    Channel
        .fromPath("${params.fastq_dir}/*_interleaved.fastq.gz")
        .map { fastq ->
            def id = fastq.name.replaceFirst(/_interleaved\.fastq\.gz$/, '')
            tuple(id, fastq)
            //tuple(fastq.baseName.replaceFirst(/_interleaved$/, ''), fastq)
        }
        .set { fastq_ch }

    bams = align_reads(fastq_ch)

    bams_ch = bams
        .map { tuple -> tuple[1] }
        .collect()
    
    merge_bam(bams_ch) 
        | variant_calling
//    variant_calling(bams_ch)
}

process align_reads {
    input:
    tuple val(sample_id), path(fastq)
	
    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    bwa mem -p ${params.ref} ${fastq} |
        samtools view -b - |
        samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}" - |
        samtools sort -o ${sample_id}.bam
    #bwa index ${sample_id}.bam
    """
}

process merge_bam {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path bams

    output:
    path "merged.bam"

    script:
    """
    echo "Starting merge"
    samtools merge merged.bam ${bams}
    """
}

process variant_calling {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path bam

    output:
    path "exome.vcf.gz"

    script:
    """
    echo "Indexing"
    bwa index ${bam}

    echo "Starting variant calling"
    bcftools mpileup -Ou -f ${params.ref} -a INFO/AD ${bam} |
        bcftools call -mv -Oz -o exome.vcf.gz
    """
}
