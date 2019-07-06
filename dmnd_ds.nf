#!/usr/bin/env  nextflow


dmnd_files = Channel.fromPath("/home/drewx/Documents/samsa2/bin/WD/*")
reference  = Channel.value("/home/drewx/Documents/samsa2/bin/extract_Combined_protozoa_bacteria.protein.faa")
OUT_DIR    = "${PWD}/dmnd.out"


process  dmnd_analysis{

    echo true
    cpus 1
    memory "2 GB"
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