#Download DNase, Benzonase, Cyanase, MNase data
fasterq-dump SRR535735
fasterq-dump SRR535736
fasterq-dump SRR535737
fasterq-dump SRR535738
fasterq-dump SRR535739
fasterq-dump SRR535740
fasterq-dump SRR535741
fasterq-dump SRR535742
fasterq-dump SRR535743
fasterq-dump SRR535744
fasterq-dump SRR5723785

#Rename DNase, Benzonase, Cyanase, MNase data from SRR number
mv SRR535735.fastq DNaseI_a.fastq
mv SRR535736.fastq DNaseI_b.fastq
mv SRR535737.fastq mm39_liver_Benzonase0.25U.fastq
mv SRR535738.fastq mm39_liver_Benzonase1U_1.fastq
mv SRR535739.fastq mm39_liver_Benzonase1U_2.fastq
mv SRR535740.fastq mm39_liver_Benzonase4U.fastq
mv SRR535741.fastq mm39_liver_Cyanase0.25U.fastq
mv SRR535742.fastq mm39_liver_Cyanase1U_1.fastq
mv SRR535743.fastq mm39_liver_Cyanase1U_2.fastq
mv SRR535744.fastq mm39_liver_Cyanase4U.fastq
mv SRR5723785_1.fastq mm39_liver_MNase_1.fastq
mv SRR5723785_2.fastq mm39_liver_MNase_2.fastq

#Concatenate each enzyme
cat DNaseI_*.fastq > mm39_liver_DNase.fastq
cat *Benz* > mm39_liver_Benzonase.fastq
cat *Cyan* > mm39_liver_Cyanase.fastq
cat *MNase* > mm39_liver_MNase.fastq

#Zip each file
gzip *ase.fastq

#build index genome for mm39
bowtie2-build mm39.fna mm39

#Align each dataset to mm39 genome
for fq in mm39_liver_*.fastq.gz
do
    name=$(echo $fq | awk -F".fastq.gz" '{print $1}')
    echo $name
    bowtie2 -p 6 --maxins 500 -x mm39 -U ${name}.fastq.gz  \
       | samtools view -bS - | samtools sort - -o $name.bam
done

#First, split each dataset into plus and minus aligned reads.
#Then, run seqOutBias on plus/minus and unseparated strands to get
##the bedfile which has chromosome coordinates for each read.
#From the bedfile, make each read a unique entry using bedToOneEntryBed.py
#Using awk commands, shift each position back 10bp then forward 21 in order to
##make a 'window' around the cut site.
#In order to get each sequence, we then convert from fasta to bed format for each data set
for i in mm39*ase.bam
do
    name=$(echo $i | awk -F".bam" '{print $1}')
    echo $name
    samtools view -bh -F 20 ${name}.bam > ${name}_plus.bam
    samtools view -bh -f 0x10 ${name}.bam > ${name}_minus.bam
    seqOutBias mm39.fna ${name}.bam --no-scale --shift-counts \
                                 --bed=${name}.bed \
                                 --bw=${name}.bigWig --read-size=60
    seqOutBias mm39.fna ${name}_plus.bam --no-scale --shift-counts \
                                 --bed=${name}_plus.bed \
                                 --bw=${name}_plus.bigWig --read-size=60
    seqOutBias mm39.fna ${name}_minus.bam --no-scale --shift-counts \
                                 --bed=${name}_minus.bed \
                                 --bw=${name}_minus.bigWig --read-size=60
    python bedToOneEntryBed.py -i ${name}.bed
    python bedToOneEntryBed.py -i ${name}_plus.bed
    python bedToOneEntryBed.py -i ${name}_minus.bed
#these bed files can be used to get the sequence flanking ALL Tn5 insertion sites,
#so we can see the precise nature of the sequence bias for each PE/strand combination
    awk '{$2 = $2 - 10; print}' ${name}.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 21; print}' | grep -v - | \
       fastaFromBed -fi mm39.fna -s -bed stdin -fo ${name}.fasta
    awk '{$2 = $2 - 10; print}' ${name}_plus.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 21; print}' | grep -v - | \
       fastaFromBed -fi mm39.fna -s -bed stdin -fo ${name}_plus.fasta
    awk '{$2 = $2 - 11; print}' ${name}_minus.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 21; print}' |  grep -v - | \
       fastaFromBed -fi mm39.fna -s -bed stdin -fo ${name}_minus.fasta
done

#First convert all sequences in the minus fasta into uppercase letters,
#then reverse complement all minus strand sequences.
for i in mm39*ase.bam
do
    name=$(echo $i | awk -F".bam" '{print $1}')
    echo $name
    awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${name}_minus.fasta | \
        fastx_reverse_complement -o ${name}_minus_RC.fasta
done
