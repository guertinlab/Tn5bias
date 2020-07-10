

#comments are for Rivanna (in progress)
wget https://raw.githubusercontent.com/guertinlab/Tn5bias/master/bedToOneEntryBed.py
#module load gcc/7.1.0
#module load seqoutbias/1.2.0
#module load mvapich2/2.3.1
#module load openmpi/3.1.4
#module load gcc/system
#module load intel/18.0  
#module load bedtools/2.26.0
#module load python/3.7.7
#ijob -A guertinlab -c 1 -p standard -t 10:00:00
for i in C1_gDNA_rep2.bam
do
    name=$(echo $i | awk -F".bam" '{print $1}')
    echo $name
#separate PE1 and PE2 reads
    samtools view -b -f 0x0040 ${name}.bam > ${name}_PE1.bam
    samtools view -b -f 0x0080 ${name}.bam > ${name}_PE2.bam
#separate plus and minus strand aligning reads:
    samtools view -bh -F 20 ${name}_PE1.bam > ${name}_PE1_plus.bam
    samtools view -bh -f 0x10 ${name}_PE1.bam > ${name}_PE1_minus.bam
    samtools view -bh -F 20 ${name}_PE2.bam > ${name}_PE2_plus.bam
    samtools view -bh -f 0x10 ${name}_PE2.bam > ${name}_PE2_minus.bam
#run seqOutBias to get the bed
    seqOutBias hg38.fa ${name}_PE1_plus.bam --no-scale --out=no_scale.tbl --shift-counts --bw=${name}_PE1_plus.bigWig --read-size=30
    seqOutBias hg38.fa ${name}_PE1_minus.bam --no-scale  --out=no_scale.tbl --shift-counts --bw=${name}_PE1_minus.bigWig --read-size=30
    seqOutBias hg38.fa ${name}_PE2_plus.bam --no-scale  --out=no_scale.tbl --shift-counts --bw=${name}_PE2_plus.bigWig --read-size=30
    seqOutBias hg38.fa ${name}_PE2_minus.bam --no-scale --out=no_scale.tbl --shift-counts --bw=${name}_PE2_minus.bigWig --read-size=30
    python bedToOneEntryBed.py ${name}_PE1_plus_not_scaled.bed
    python bedToOneEntryBed.py ${name}_PE1_minus_not_scaled.bed
    python bedToOneEntryBed.py ${name}_PE2_plus_not_scaled.bed
    python bedToOneEntryBed.py ${name}_PE2_minus_not_scaled.bed
#these bed files can be used to get the sequence flanking ALL Tn5 insertion sites, 
#so we can see the precise nature of the sequence bias for each PE/strand combination
    awk '{$2 = $2 - 9; print}' ${name}_PE1_plus_not_scaled.bed | awk '{OFS="\t";} {$3 = $2 + 20; print}' | grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_PE1_plus_not_scaled.fasta
    awk '{$2 = $2 - 9; print}' ${name}_PE1_minus_not_scaled.bed | awk '{OFS="\t";} {$3 = $2 + 20; print}' |  grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_PE1_minus_not_scaled.fasta
    awk '{$2 = $2 - 9; print}' ${name}_PE2_plus_not_scaled.bed | awk '{OFS="\t";} {$3 = $2 + 20; print}' | grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_PE2_plus_not_scaled.fasta
    awk '{$2 = $2 - 9; print}' ${name}_PE2_minus_not_scaled.bed | awk '{OFS="\t";} {$3 = $2 + 20; print}' |  grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_PE2_minus_not_scaled.fasta
done
