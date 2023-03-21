# seqOutATACBias: Rule Ensemble Modeling of seqOutBias Scaling
[![DOI](https://zenodo.org/badge/617080961.svg)](https://zenodo.org/badge/latestdoi/617080961)    
This GitHub repo contains all of the code and information produced for the Guertinlab Tn5 bias correction project. It was made to accompany the paper "Correction of transposase sequence bias in ATAC-seq data with rule ensemble modeling" by providing the code and methods used to produce the figures and data. It is divided into 4 sections:


## seqOutATACBias_setup
The folder labeled 'seqOutATACBias_setup' contains seqOutATACBias, which is a CLI that corrects the sequence bias of Tn5 transposase in ATAC-seq data using a rule ensemble model.


## seqOutATACBias_workflow_Vignette
The folder labeled 'seqOutATACBias_workflow_Vignette' contains a vignette which explains the methods used by seqOutATACBias and uses data from chromosome 21 to show bias correction. This is a light weight analysis that can be conducted in less than 15 minutes, and uses less than 1Gb of disk space when completed.    


## Manuscript_Vignette
The folder labled 'Manuscript_Vignette' contains all of the code to conduct the analysis outlined in "Correction of transposase sequence bias in ATAC-seq data with rule ensemble modeling". This vignette will download all of the publicly available data used in the paper and conduct the appropriate analysis. Completion takes several days and requires the use of an HPC environment and slightly more than a terabyte of storage space.


## Manuscript_Figures
The 'Manuscript_Figures' folder contains the data which is directly plotted in "Correction of transposase sequence bias in ATAC-seq data with rule ensemble modeling" figures and the necessary code to produce said figures.
