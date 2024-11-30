# 5hmC
The codes shown here were used to perform the alignment of the fastq files to human genome (hg38), peak calling, elastic net selection, and regression for identifying AD related 5hmCs.

## analysis process
Alignment:

    sequencing_process.sh


Peak calling:

    peakCalling_paired.sh
    peakCalling_single.sh
    derive consensus peaks.R


Elastic Net selection:

    elastic net model.R
    variation explained.R

    
Regression:

    association analysis.


Colocalization analysis:

    Colocalization Analysis.R


R code for the figure of demographics of multi-omics data:

    Figure--Demographics of Multi-omics Data.R


GARFIELD shell script:

    GARFIELD_Code.sh

