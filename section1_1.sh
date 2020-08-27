fasterq-dump SRR5123141
fasterq-dump SRR5123142

gzip *fastq

mv SRR5123141_1.fastq.gz C1_gDNA_rep1_PE1.fastq.gz
mv SRR5123141_2.fastq.gz C1_gDNA_rep1_PE2.fastq.gz
mv SRR5123142_1.fastq.gz C1_gDNA_rep2_PE1.fastq.gz
mv SRR5123142_2.fastq.gz C1_gDNA_rep2_PE2.fastq.gz

#you also need to download the human genome:
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chrM.fa.gz

#unzip
gunzip hg38.fa.gz
gunzip chrM.fa.gz
