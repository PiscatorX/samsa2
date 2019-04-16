#!/usr/bin/env nextflow

SAMSA			= "/home/drewx/Documents/samsa2/"
params.INPUT_DIR	= "/home/drewx/Documents/samsa2/sample_files_paired-end/1_starting_files/*_R{1,2}.fastq"
params.OUT_DIR		= "${PWD}/samsa2Out"
OUT_DIR                 = params.OUT_DIR
params.diamond_refseq	= "$SAMSA/setup_and_test/tiny_databases/RefSeq_bac_TINY_24MB"
params.diamond_subsys_db= "$SAMSA/setup_and_test/tiny_databases/subsys_db_TINY_24MB"
params.RefSeq_db	= "$SAMSA/setup_and_test/tiny_databases/RefSeq_bac_TINY_24MB.fa"
params.Subsys_db	= "$SAMSA/setup_and_test/tiny_databases/subsys_db_TINY_24MB.fa"
params.sortmerna_fasta 	= "/home/drewx/Documents/samsa2/programs/sortmerna-2.1/rRNA_databases/silva-bac-16s-id90.fasta"
params.sortmerna_index  = "/home/drewx/Documents/samsa2/programs/sortmerna-2.1/index/silva-bac-16s-db"
sortmerna_fasta         =  Channel.value(params.sortmerna_fasta)
sortmerna_index         =  Channel.value(params.sortmerna_index)
diamond_refseq          =  Channel.value(params.diamond_refseq)
diamond_subsys_db       =  Channel.value(params.diamond_subsys_db) 
RefSeq_db	        =  Channel.value(params.RefSeq_db)	
Subsys_db	        =  Channel.value(params.Subsys_db)	


Channel.fromFilePairs(params.INPUT_DIR)
           .ifEmpty{ error "Could not locate pair reads: ${reads}"}
           .into{reads1; reads2}



process trimmomatic{

    //echo true
    cpus params.htp_cores
    memory "${params.l_mem} GB"
    publishDir path: "${OUT_DIR}/trimmomatic", mode: 'copy'
    
    input:
	set pair_id, file(reads) from reads1  

    output:
	set pair_id, file("${pair_id}.forward"), file("${pair_id}.reverse") into  trimmmomatic_reads1
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
    ${pair_id}.forward \
    ${pair_id}.forward_unpaired \
    ${pair_id}.reverse \
    ${pair_id}.reverse_unpaired \
    SLIDINGWINDOW:4:15 MINLEN:70

"""
    
}



process raw_read_counter{

    memory "${params.l_mem} GB"
    //echo true    
    input:
        set pair_id, file(trimm_fwd), file(trimm_rev) from trimmmomatic_reads1
    output:
	 set pair_id, file(trimm_fwd), file(trimm_rev) into trimmmomatic_reads2


"""    
 
  raw_read_counter.py \
  -I ${trimm_fwd} \
  -O ${OUT_DIR}/trimmomatic/raw_counts.txt
 
"""

}




process pear{
	
    cpus params.ltp_cores
    memory "${params.l_mem} GB"
    publishDir path: "${OUT_DIR}/pear/pair_id", mode: 'copy'

    input:
        set pair_id, file(trimm_fwd), file(trimm_rev) from trimmmomatic_reads2

    output:
         set pair_id, file("${pair_id}.assembled.fastq") into merged_reads
         set file("${pair_id}.discarded.fastq"),  file("${pair_id}.unassembled*") into umerged_reads
    
"""

   pear \
   -f ${trimm_fwd} \
   -r ${trimm_rev} \
   --memory ${params.l_mem}G \
   --threads ${params.htp_cores} \
   -o ${pair_id}
  
"""


}




process  sortmerna{

    echo true
    cpus params.ltp_cores
    memory "${params.l_mem} GB"
    publishDir path: "${OUT_DIR}/sortmerna", mode: 'copy'
    
    input:
        set pair_id, file(assembled_reads) from merged_reads       
	val sortmerna_fasta
	val sortmerna_index

    output:
	file("*.ribosomes.fastq*") into ribosomes
        file("${pair_id}.ribosomes.log") into  sortmerna_log
        file("${pair_id}.fastq") into metatranscriptome_reads
        val  pair_id   

"""

    sortmerna \
    -a ${params.htp_cores} \
    --ref ${sortmerna_fasta},${sortmerna_index} \
    --reads ${assembled_reads} \
    --aligned ${pair_id}.ribosomes \
    --other ${pair_id} \
    --fastx \
    --log 
  
"""
    
}




process diamond{

    //echo true
    cpus params.ltp_cores
    memory "${params.m_mem} GB"
    publishDir path: "${OUT_DIR}/diamond", mode: 'copy'
    input:
	file(query_seqs) from metatranscriptome_reads
        val  diamond_refseq
        val  pair_id


   output:
        file("${pair_id}.daa") into refseq_daa
	file("${pair_id}.tab") into refseq_tab
    

"""
 
    diamond \
    blastx \
    -d ${diamond_refseq}  \
    -q ${query_seqs} \
    -a ${pair_id}.daa \
    --max-target-seqs 1 \
    --verbose

    diamond \
    view \
    --daa ${pair_id}.daa \
    -o ${pair_id}.tab \
    -f tab    

"""
    
   
}


process analysis_counter{

    echo true
    input:
       file refseq_results from refseq_tab




"""
   
   DIAMOND_analysis_counter.py -I ${refseq_results} -D $RefSeq_db -O
   DIAMOND_analysis_counter.py -I ${refseq_results} -D $RefSeq_db -F
   tree
 
"""

}































//daa = Channel.fromPath("/home/drewx/Documents/SAMSA-TEST/samsa2/output_test/step_4_output_test/daa_binary_files/*.daa")





// process diamond{

//     maxForks 1
//     echo true
//     input:

//         file(daa_file) from daa
// 	//file(query_seqs) from metatranscriptome_reads.flatten()
//         //val  diamond_db
//         //val  pair_id
    

// """

   
//     diamond \
//     view \
//     --daa ${daa_file} \
//     -o ${daa_file}.tab \
//     -f tab    
 
// """
    
   
// }



// echo diamond \
// blastx \
// -d ${diamond_db}  \
// -q ${query_seqs} \
// -a ${pair_id}.daa \
// --block-size 0.1 \
// --index-chunks 20 \
// --max-target-seqs 1 \
// --verbose

// diamond \
// blastx \
// -d ${diamond_db}  \
// --un ${pair_id}.unaligned \
// --al ${pair_id}.aligned \
// -q ${query_seqs} \
// -a ${pair_id}.daa \
// -f ${params.outformat}  \
// --more-sensitive \
// --frameshift 15 \
// --max-target-seqs 1 \
// --verbose    
