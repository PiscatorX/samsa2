#!/usr/bin/env nextflow

SAMSA="/home/drewx/Documents/samsa2/"
params.INPUT_DIR="/home/drewx/Documents/samsa2/sample_files_paired-end/1_starting_files/*_R{1,2}.fastq"
params.OUT_DIR="${PWD}/samsa2Out"
//params.getconf _NPROCESSORS_ONLN= getconf _NPROCESSORS_ONLN
OUT_DIR=params.OUT_DIR
params.diamond_database="$SAMSA/setup_and_test/tiny_databases/RefSeq_bac_TINY_24MB"
params.diamond_subsys_db="$SAMSA/setup_and_test/tiny_databases/subsys_db_TINY_24MB"
params.RefSeq_db="$SAMSA/setup_and_test/tiny_databases/RefSeq_bac_TINY_24MB.fa"
params.Subsys_db="$SAMSA/setup_and_test/tiny_databases/subsys_db_TINY_24MB.fa"
STEP_1="$OUT_DIR}/step_1_output"
STEP_2="$OUT_DIR}/step_2_output"
STEP_3="$OUT_DIR}/step_3_output"
STEP_4="$OUT_DIR}/step_4_output"
STEP_5="$OUT_DIR}/step_5_output"




Channel.fromFilePairs(params.INPUT_DIR)
           .ifEmpty{ error "Could not locate pair reads: ${reads}"}
           .into{reads1; reads2}



process trimmomatic{

    echo true
    cpus params.htp_cores
    memory "${params.h_mem} GB"
    publishDir path: "${OUT_DIR}/trimmomatic", mode: 'copy'
    
    input:
	set out_path, file(reads) from reads1  

    output:
	set file("${out_path}.forward"), file("${out_path}.reverse") into  (trimmmomatic_reads1, trimmmomatic_reads2)    
        file("*_unpaired") into trimmed__unpaired     

    script:
	(fwd, rev) = reads
    
    
"""

    $trimmomatic \
    PE \
    -phred33 \
    -threads ${params.htp_cores} \
    ${fwd} \
    ${rev} \
    ${out_path}.forward \
    ${out_path}.forward_unpaired \
    ${out_path}.reverse \
    ${out_path}.reverse_unpaired \
    SLIDINGWINDOW:4:15 MINLEN:70

"""
    
}



process raw_read_counter{

    maxForks 1
    echo true    
    input:
        set file(trimm_fwd), file(trimm_rev) from trimmmomatic_reads1



"""    
 
  raw_read_counter.py \
  -I ${trimm_fwd} \
  -O ${OUT_DIR}/trimmomatic/raw_counts.txt
 
"""

}




// process pear{
	
//     //echo true
//     cpus params.htp_cores
//     memory "${params.h_mem} GB"
//     publishDir path: OUT_DIR, mode: 'copy'

//     input:
//         set out_path, file(trimm_fwd), file(trimm_rev) from trimmmomatic_reads2

    
// """

//    pear \
//    -f ${trimm_fwd} \
//    -r ${trimm_rev} \
//    --memory ${params.h_mem}G \
//    --threads ${params.htp_cores} \
//    -o out_path

// """							    							 
// }


