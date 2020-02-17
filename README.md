# SAMSA2-nf: a nextflow implementation of the samsa2 metatranscriptome analysis pipeline

*****
For more information about  samsa2 please visit the repository (<https://github.com/transcript/samsa2>)  with the original samsa2 implementation.

The approach of the SAMSA2 pipeline was followed in the implementation of the  pipeline. The pipeline has been implemented in nexflow instead of shell scripting. Not all functions have been implicated.

## Prerequisites

### Nextflow

Nextflow requires Java to to run more details from [https://www.nextflow.io/](https://www.nextflow.io/)

### Dependincies

Dependencies are listed in the original samsa2 repository [here](https://github.com/transcript/samsa2#dependencies).
Add dependencies to the the PATH environmental variable.

### Running the pipeline

The main nextflow script is [samsa2.nf](https://github.com/PiscatorX/samsa2-nf/blob/master/samsa2.nf) and can be run by the command ``nextflow  samsa2.nf``

### References

Westreich, S. T., Treiber, M. L., Mills, D. A., Korf, I., & Lemay, D. G. (2018). SAMSA2: a standalone metatranscriptome analysis pipeline. BMC bioinformatics, 19(1), 175. [https://dx.doi.org/10.1186%2Fs12859-018-2189-z](https://dx.doi.org/10.1186%2Fs12859-018-2189-z) 