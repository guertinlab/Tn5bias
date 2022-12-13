# seqOutATACBias Workflow Vignette

In this streamlined vignette, an example rule ensemble correction is applied in order to display the workflow used in **Correction of transposase sequence bias in ATAC-seq data with rule ensemble modeling**.
Using a subset of the data in the paper, reads aligning to and the reference genome sequences for chromosome 21 are first downloaded. 
seqOutBias is ran to generate inputs for the rule ensemble model. 
A pretrained rule ensemble model is then applied to these inputs. 
Finally, the rule ensemble and seqOutBias output is plotted using the coordinates for ESR1 on chromosome 21.

Please check that your local system has all necessary dependencies before running the RMD file.

**Estimated run time: 15 minutes**       
**Size on disk: less than 650Mb**      
