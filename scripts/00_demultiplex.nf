nextflow.enable.dsl=2

// baseline parameters
params.basedir = "/90daydata/guedira_seq_map"
params.fastq_dir = "${params.basedir}/fastq"
params.keyfile = "/project/guedira_seq_map/nico/SunRILs_population_description/data/keyfiles/SNP_calling_cladoclean_aviti_illumina_20250817.tsv"
//params.keyfile = "${params.basedir}/nico2/test_keyfile.tsv"
params.output_dir = "${params.basedir}/nico2/demultiplexed_fastq"


workflow {
    // make output directory if it doesn't exist
    new File(params.output_dir).mkdirs()


    Channel
        .fromPath(params.keyfile)
	.splitCsv(header: true, sep: '\t')
            .map { row ->
	        tuple(
		    row.FullSampleName,
		    row.Barcode,
		    row.Flowcell,
		    row.Lane
		)
	    }
	.set { key_ch }
	
    Channel
        .fromPath("${params.fastq_dir}/*fastq.gz")
	.map { fq ->
            def name = fq.baseName
            def parts = name.tokenize('_')

            assert parts.size() >= 2 : "Unexpected FASTQ name: ${name}"

            def flowcell = parts[0]
            def lane = parts[1]

	    tuple(flowcell, lane, fq)
	}
	.set { fastq_ch }

    key_ch
        .join(fastq_ch, by: [2,3])
	.map { sample, barcode, flowcell, lane, fq ->
		tuple(sample, barcode, fq)
	}
	.set { demux_jobs }

    key_ch.view { "KEY → $it" }
    demux_jobs.view { "JOINED → $it" }

    grouped_samples = (
	    key_ch 
	    | CUTADAPT_DEMULTIPLEX
	    | groupTuple(by: 0)
	)
	
    grouped_samples | MERGE_SAMPLE_FASTQS
}


process CUTADAPT_DEMULTIPLEX {
    input:
    tuple val(sample), val(barcode), val(flowcell), val(lane)
   
    output:
    tuple val(sample), path("${sample}_${flowcell}_${lane}_${barcode}.fastq.gz")	
    
    script:
    """
    cutadapt \
        --cores=0 \
        -g ^${barcode} \
        --discard-untrimmed \
        -e 0 \
        --quality-cutoff 30 \
        --minimum-length 50 \
        --max-n 0 \
        -o ${sample}_${flowcell}_${lane}_${barcode}.fastq.gz \
        ${params.fastq_dir}/${flowcell}_${lane}_fastq.gz
    """
}

process MERGE_SAMPLE_FASTQS {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple val(sample), path(fqs)

    output:
    path "${sample}.fastq.gz"

    script:
    """
    cat ${fqs.join(' ')} > ${sample}.fastq.gz
    zcat *.fastq.gz | echo \$((`wc -l` / 4))
    """
}
