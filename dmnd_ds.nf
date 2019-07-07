#!/usr/bin/env  nextflow


dmnd_files = Channel.fromPath("/home/andhlovu/MT-samsa2/samsa2.Out/diamond_refseq/*.tab")
reference  = Channel.value("/projects/andhlovu/DB_REF/RefSeq/Combined_protozoa_bacteria.protein.faa")
OUT_DIR    = "${PWD}/dmnd.out"


process  dmnd_analysis{

    echo true
    cpus 2
    memory "10 GB"
    publishDir path: OUT_DIR, mode: 'move'
    input:
       val reference
       file dmnd from dmnd_files


    output:
       file("*")
    
    script:
       prefix =  "SHB_"


"""    

   dmnd_analysis.py \
   ${dmnd} \
   -r ${reference}\
   -p ${dmnd}

"""


}