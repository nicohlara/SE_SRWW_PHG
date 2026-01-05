nextflow.enable.dsl=2

//to run must have nextflow and miniconda3 loaded.

//  Set up parameters
params.input_file = file('/90daydata/guedira_seq_map/nico/plinkhaplo_phg/output/slurm_align_file.txt')

Channel
    .fromPath(params.input_file)
    .splitText()
    .set { command_lines }

process phg_align {
    conda '/home/nicolas.lara/.conda/envs/phgv2-conda'
   // publishDir params.output_dir, mode: 'copy'

    input:
    val cmd_line
	
    script:
    """
    cd /90daydata/guedira_seq_map/nico/pangenome_multichrom/
    echo "Running: ${cmd_line}"
    eval ${cmd_line} 
    """
}

workflow {
    phg_align(command_lines)
}

