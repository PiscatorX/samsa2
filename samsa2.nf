#!/bin/env nextflow

SAMSA="/home/drewx/Documents/samsa2/"
params.INPUT_DIR=/home/drewx/Documents/samsa2/sample_files_paired-end/1_starting_files/
params.OUT_DIR="${PWD}/sams2Out"
params.getconf _NPROCESSORS_ONLN= `getconf _NPROCESSORS_ONLN`
params.outdir=OUT_DIR
params.diamond_database="$SAMSA/setup_and_test/tiny_databases/RefSeq_bac_TINY_24MB"
params.diamond_subsys_db="$SAMSA/setup_and_test/tiny_databases/subsys_db_TINY_24MB"
//ggregation databases
params.RefSeq_db="$SAMSA/setup_and_test/tiny_databases/RefSeq_bac_TINY_24MB.fa"
params.Subsys_db="$SAMSA/setup_and_test/tiny_databases/subsys_db_TINY_24MB.fa"
STEP_1="$OUT_DIR}/step_1_output"
STEP_2="$OUT_DIR}/step_2_output"
STEP_3="$OUT_DIR}/step_3_output"
STEP_4="$OUT_DIR}/step_4_output"
STEP_5="$OUT_DIR}/step_5_output"



process Step1_trimmomatic{

"""

    java -jar $TRIMMOMATIC \ 
    PE \
    -phred33 \
    -threads $threads \
    $f \
    $f2  
    $out_path".forward" \
    $out_path".forward_unpaired" \
    $out_path".reverse" \
    $out_path".reverse_unpaired" \
    SLIDINGWINDOW:4:15 MINLEN:70

"""

}
