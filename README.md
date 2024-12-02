# 5hmC associated with Alzheimer's disease (Brain 5-hydroxymethylcytosine alterations are associated with Alzheimer’s disease neuropathology)
The programs shown here were used to perform the alignment of the fastq files to human genome (hg38), peak calling, and downstream analyses (e.g., association analysis, enrichment analysis) for identifying AD neuropathology-associated 5hmCs.

## Statistical analysis
5hmC pre-processing:

    sequencing_process.sh
    peakCalling_paired.sh
    peakCalling_single.sh
    derive consensus peaks.R


Elastic Net selection:

    elastic net model.R

    
Association analysis:

    association analysis.R
    variation explained.R
    meta-analysis.R
    

Genomic feature analysis:

    genomic feature plot.R


Colocalization analysis:

    Colocalization Analysis.R


R code for the figure of demographics of multi-omics data:

    Figure--Demographics of Multi-omics Data.R


GARFIELD shell script:

    GARFIELD_Code.sh

