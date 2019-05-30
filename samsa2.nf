#!/usr/bin/env nextflow

//params.INPUT_DIR        = "/home/andhlovu/Novogene/ftpdata.novogene.cn:2300/C101HW18111065/raw_data/*_RNA_{1,2}.fq.gz" 
params.INPUT_DIR        =  "/home/drewx/Documents/samsa2/samsa2.Out/sortmerna/*.fastq"


params.OUT_DIR          = "${PWD}/samsa2.Out"
OUT_DIR                 = params.OUT_DIR

// params.diamond_subsys_db= "/projects/andhlovu/DB_REF/Subsys/subsys_db.dmnd"
// params.Subsys_db        = "/projects/andhlovu/DB_REF/Subsys/subsys_db.fa"
// params.diamond_refseq   = "/projects/andhlovu/DB_REF/RefSeq/RefSeq_bac.dmnd"
// params.RefSeq_db        = "/projects/andhlovu/DB_REF/RefSeq/RefSeq_bac.fa"
//params.sortmerna_fasta  = "/projects/andhlovu/DB_REF/SILVA/SILVA_132_SSURef_Nr99_tax_silva.fasta"
//params.sortmerna_index  = "/projects/andhlovu/DB_REF/SortMeRNA/SILVA.idx"

params.diamond_refseq	   = "/home/drewx/Documents/SAMSA-TEST/samsa2/setup_and_test/tiny_databases/RefSeq_bac_TINY_24MB"
params.diamond_subsys_db   = "/home/drewx/Documents/SAMSA-TEST/samsa2/setup_and_test/tiny_databases/subsys_db_TINY_24MB"
params.RefSeq_db 	   = "/home/drewx/Documents/SAMSA-TEST/samsa2/setup_and_test/tiny_databases/RefSeq_bac_TINY_24MB.fa"
params.Subsys_db 	   = "/home/drewx/Documents/SAMSA-TEST/samsa2/setup_and_test/tiny_databases/subsys_db_TINY_24MB.fa"
params.sortmerna_fasta     = "/home/drewx/Documents/samsa2/programs/sortmerna-2.1/rRNA_databases/silva-bac-16s-id90.fasta"
params.sortmerna_index     = "/home/drewx/Documents/samsa2/programs/sortmerna-2.1/index/silva-bac-16s-db"


sortmerna_fasta         =  Channel.value(params.sortmerna_fasta)
sortmerna_index         =  Channel.value(params.sortmerna_index)
diamond_refseq          =  Channel.value(params.diamond_refseq)
diamond_subsys_db       =  Channel.value(params.diamond_subsys_db)
RefSeq_db               =  Channel.value(params.RefSeq_db)
Subsys_db               =  Channel.value(params.Subsys_db)
params.diamond_only     =  true


if(! params.diamond_only ){

Channel.fromFilePairs(params.INPUT_DIR)
           .ifEmpty{ error "Could not locate pair reads: ${params.INPUT_DIR}"}
           .into{reads1; reads2}


process trimmomatic{

    //echo true
    cpus params.mtp_cores
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

    maxForks 1
    cpus params.ltp_cores
    memory "${params.l_mem} GB"
    //echo true    
    input:
        set pair_id, file(trimm_fwd), file(trimm_rev) from trimmmomatic_reads1
    output:
	 set pair_id, file(trimm_fwd), file(trimm_rev) into trimmmomatic_reads2


"""    
 
  raw_read_counter.py \
  -I ${trimm_fwd} \
  -O ${OUT_DIR}/raw_counts.txt
 
"""

}




process pear{
	
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
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

    //echo true
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
    publishDir path: "${OUT_DIR}/sortmerna", mode: 'copy'
    
    input:
        set pair_id, file(assembled_reads) from merged_reads       
	val sortmerna_fasta
	val sortmerna_index

    output:
	file("*.ribosomes.fastq*") into ribosomes
        file("${pair_id}.ribosomes.log") into  sortmerna_log
        file("${pair_id}.fastq") into (metatranscriptome_reads1, metatranscriptome_reads2)
        val  pair_id into (pair_id1, pair_id2)  

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

}

if (params.diamond_only ){

Channel.fromPath(params.INPUT_DIR)
           .ifEmpty{ error "Could not metranscriptome reads : ${params.INPUT_DIR}"}
           .into{get_pair_id; metatranscriptome_reads1; metatranscriptome_reads2}
	   
get_pair_id.map{it.baseName}.into{pair_id1; pair_id2}

}


process diamond_Refseq{

    //echo true
    cpus params.htp_cores
    memory "${params.m_mem} GB"
    publishDir path: "${OUT_DIR}/diamond_refseq", mode: 'copy'
    input:
	file(query_seqs) from metatranscriptome_reads1
        val  diamond_refseq
        val  pair_id from pair_id1


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


process diamond_subsys{

    //echo true
    cpus params.htp_cores
    memory "${params.m_mem} GB"
    publishDir path: "${OUT_DIR}/diamond_subsys", mode: 'copy'
    input:
	file(query_seqs) from metatranscriptome_reads2
        val diamond_subsys_db  
        val  pair_id from pair_id2

   output:
        file("${pair_id}.daa") into subsys_daa
	file("${pair_id}.tab") into subsys_tab
    

"""
 
    diamond \
    blastx \
    -d ${diamond_subsys_db}  \
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

    //echo true
    publishDir path: "${OUT_DIR}/RefSeq_results", mode: 'copy'
    input:
       file refseq_results from refseq_tab
       val RefSeq_db   

    output:
        file("*organism.tsv") into refseq_org_results
        file("*function.tsv") into refseq_func_results

    script:
	outfile = refseq_results.getName().replaceFirst(/tab/, "stdout") 
    

"""
  
   mkdir -pv ${OUT_DIR}/RefSeq_top_results
   DIAMOND_analysis_counter.py -I ${refseq_results} -D ${RefSeq_db} -O > ${OUT_DIR}/RefSeq_top_results/orgnism_${outfile}  
   DIAMOND_analysis_counter.py -I ${refseq_results} -D ${RefSeq_db} -F > ${OUT_DIR}/RefSeq_top_results/function_${outfile}
 
"""

}


process analysis_counter_subsys{

    //echo true
    cpus params.ltp_cores
    memory "${params.l_mem} GB"
    publishDir path: "${OUT_DIR}/Subsys_results", mode: 'copy'
    input:
       file subsys_results from subsys_tab
       val Subsys_db   

    output:
	 file("*.reduced")  into sub_sys_reduced
	 file("*.receipt") into sub_sys_receipt
         file("*.hierarchy") into sub_sys_hierarchy
    script:
         outfile = subsys_results.getName().replaceFirst(/tab/, "stdout") 
    
    
    
"""
   
   mkdir -pv ${OUT_DIR}/Subsys_top_results
   
   DIAMOND_subsystems_analysis_counter.py \
   -I ${subsys_results} \
   -D ${Subsys_db} \
   -O ${subsys_results}.hierarchy \
   -P ${subsys_results}.receipt > ${OUT_DIR}/Subsys_top_results/${outfile}

   subsys_reducer.py \
  -I ${subsys_results}.hierarchy   

 
"""

}
