# 5hmC associated with Alzheimer's disease (Brain 5-hydroxymethylcytosine alterations are associated with Alzheimer’s disease neuropathology)
The programs shown here were used to perform the alignment of the fastq files to human genome (hg38), peak calling, and downstream analyses (e.g., association analysis, enrichment analysis) for identifying AD neuropathology-associated 5hmCs.

## Statistical analysis
5hmC pre-processing:

    sequencing_process.sh
    peakCalling_paired.sh
    peakCalling_single.sh
    generation of consensus 5hmC region.R
    batch removal.R
    PCA.R

    
Association analysis:

    elastic net model.R
    association analysis.R
    variation explained.R
    meta-analysis.R
    

Genomic feature analysis:

    genomic feature plot.R


Multi-omics integration:

    RNA-seq gene expression.R 
    ChIP-seq histone.R
    TMT proteomics.R
    

Colocalization analysis:

    Colocalization Analysis.R


GARFIELD shell script:

    GARFIELD_Code.sh

