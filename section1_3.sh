module load sratoolkit/2.9.1
fasterq-dump SRR535737
fasterq-dump SRR535738
fasterq-dump SRR535739
fasterq-dump SRR535740
fasterq-dump SRR535741
fasterq-dump SRR535742
fasterq-dump SRR535743
fasterq-dump SRR535744
fasterq-dump SRR535735
fasterq-dump SRR535736


mv SRR535737.fastq mm10_liver_Benzonase0.25U.fastq
mv SRR535738.fastq mm10_liver_Benzonase1U_1.fastq
mv SRR535739.fastq mm10_liver_Benzonase1U_2.fastq
mv SRR535740.fastq mm10_liver_Benzonase4U.fastq
mv SRR535741.fastq mm10_liver_Cyanase0.25U.fastq
mv SRR535742.fastq mm10_liver_Cyanase1U_1.fastq
mv SRR535743.fastq mm10_liver_Cyanase1U_2.fastq
mv SRR535744.fastq mm10_liver_Cyanase4U.fastq
mv SRR535735.fastq DNaseI_a.fastq 
mv SRR535736.fastq DNaseI_b.fastq

cat *Benz* > mm10_liver_Benzonase.fastq 
cat *Cyan* > mm10_liver_Cyanase.fastq
cat DNaseI_*.fastq > mm10_liver_DNase.fastq 
gzip *ase.fastq


#build index genome
#bowtie2-build mm10.fa mm10
module load gcc/7.1.0 bowtie2/2.1.0
module load samtools/1.10

for fq in mm10_liver_D*ase.fastq.gz
do
    name=$(echo $fq | awk -F".fastq.gz" '{print $1}')
    echo $name
    bowtie2 -p 4 --maxins 500 -x mm10 -U ${name}.fastq.gz  \
       | samtools view -bS - | samtools sort - -o $name.bam
#    seqOutBias mm10.fa $name.bam --no-scale --bw=${name}.bigWig --shift-counts --read-size=30
done


module load gcc/7.1.0
module load seqoutbias/1.2.0
module load fastx-toolkit/0.0.14
module load bedtools/2.26.0
module load intel/20.0  
module load intelmpi/20.0
module load python/3.7.7
#ijob -A guertinlab -c 1 -p standard -t 10:00:00

for i in mm10*ase.bam
do
    name=$(echo $i | awk -F".bam" '{print $1}')
    echo $name
#separate plus and minus strand aligning reads:
    samtools view -bh -F 20 ${name}.bam > ${name}_plus.bam
    samtools view -bh -f 0x10 ${name}.bam > ${name}_minus.bam
#run seqOutBias to get the bed
    seqOutBias mm10.fa ${name}_plus.bam --no-scale --shift-counts \
                                 --bed=${name}_plus.bed \
                                 --bw=${name}_plus.bigWig --read-size=30
    seqOutBias mm10.fa ${name}_minus.bam --no-scale --shift-counts \
                                 --bed=${name}_minus.bed \
                                 --bw=${name}_minus.bigWig --read-size=30
    python bedToOneEntryBed.py -i ${name}_plus.bed
    python bedToOneEntryBed.py -i ${name}_minus.bed
#these bed files can be used to get the sequence flanking ALL Tn5 insertion sites, 
#so we can see the precise nature of the sequence bias for each PE/strand combination
    awk '{$2 = $2 - 10; print}' ${name}_plus.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 21; print}' | grep -v - | \
       fastaFromBed -fi mm10.fa -s -bed stdin -fo ${name}_plus.fasta
    awk '{$2 = $2 - 10; print}' ${name}_minus.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 21; print}' |  grep -v - | \
       fastaFromBed -fi mm10.fa -s -bed stdin -fo ${name}_minus.fasta
done

module load gcc/9.2.0
module load fastx-toolkit/0.0.14


for i in mm10*ase.bam
do
    name=$(echo $i | awk -F".bam" '{print $1}')
    echo $name
    awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${name}_minus.fasta | \
        fastx_reverse_complement -o ${name}_minus_RC.fasta
done
