# seqOutATACBias: Rule Ensemble Modeling of seqOutBias Scaling

seqOutATACBias is a CLI which corrects the sequence bias of Tn5 transposase in ATAC-seq data.   

seqOutATACBias uses the same input as seqOutBias: a reference genome and aligned read files in BAM format.   
Output from seqOutATACBias is scaled read values in bigWig and bedGraph format.

## Current Software Verion

1.0

## Installation and setup

seqOutATACBias uses several dependencies when installing and when processing data. Please ensure that the following
software is installed on your machine and in the PATH variable:

seqOutBias >= version 1.3   
anaconda Command line client >= version 1.9.0    
python >= version 3.9.12   
setuptools python package >= version: 65.5.1   
R >= version 4.2.1   
R data.table package >= version 1.14.2   
faidx >= version 0.7.1   
GNU parallel >= 20220722   
bedtools >= version 2.30.0   
bigWigToBedGraph >= version 438   
bedGraphToBigWig >= version 2.9   

Once the required dependencies are installed and in PATH, the seqOutATACBias software can be downloaded from GitHub by entering:
```sh
$ wget https://github.com/guertinlab/Tn5bias/archive/refs/heads/master.zip
```

Upon downloading the software, you may install using either setuptools or PIP.
To install via setuptools, navigate to the seqOutATACBias_setup folder and enter:
```sh
$ python setup.py install
```

To install with PIP, first navigate to the seqOutATACBias_setup folder, then enter:
```sh
$ python setup.py sdist
$ cd dist
$ pip install seqOutATACBias-1.0.tar.gz
```

Either option will add the seqOutATACBias software to your PATH. To test that the software was correctly installed
and added to your PATH variable, enter:
```sh
$ seqOutATACBias
```

Output from a successful installation should show all declared options (or their defaults) and ask for a valid command:
```
Cleanup (-c or --cleanup)= TRUE
Read Length (-r or --readlength)= 60
Input (-i or --input)=
Genome (-g or --genome)=
Processors (-p or --processors)=
Command =
Outfile =
Please enter a command
```

## Scaling ATAC-seq data

seqOutATACBias scales ATAC-seq data by generating 12 direct k-mer scaling files and one unscaled file using seqOutBias
and applying a rule ensemble model to these files. Because seqOutATACBias requires 13 runs of seqOutBias in addition to
suffix and tallymer library construction, it is recommended that 150Gb of free disk space is alloted for each full run
(invoking the 'masks' and 'scale' commands). Each run will take approximately 12 hours of compute time, but this can be
greatly reduced by having suffix and tallymer libraries for the input reference genome and read length (seqOutATACBias will
reuse these files from previous runs or other runs of seqOutBias using the same parameters). Additionally, seqOutATACBias
will run some processes in parallel (according to the -p option), which will also reduce compute time.

To generate the input for seqOutATACBias, one must invoke the "mask" command and specify the BAM file and reference genome:
```sh
seqOutATACBias masks -i=random_scaling_test.bam -g=random.fa
```

The "Processors", "Cleanup", and "Read Length" options may also be invoked to customize the data generation process:
```sh
seqOutATACBias masks -i=random_scaling_test.bam -g=random.fa -p=5 -r=50 -c=NO
```
    
A brief description of each option:    
<br /> 
Processors: This option controls the number of processes when running seqOutBias and converting bigWigs to bedGraphs. Please ensure
that your machine has the requisite RAM and processors for your selection.    
<br /> 
Cleanup: Whether or not the .tbl, .bigWig, and .bedGraph files are deleted after their use in the next processing step.
Anything other than 'TRUE' will cause these files to not be deleted. The default is set to TRUE.
<br />    
Read Length: The length of reads in input data. Necessary for computing mappability. Default is set to 60.
<br />    
<br />    
<br />    
Once the "masks" command has finished running, an output "union" bedGraph file will be in the directory.

To apply the seqOutATACBias rule ensemble correction, the "scale" command is used with the union bedGraph file as input:
```sh
seqOutATACBias scale -i=random_scaling_test_union.bedGraph -g=random.fa
```

This process is much quicker and should only take about an hour.
Output data will be the scaled reads in bedGraph and bigWig format. 
