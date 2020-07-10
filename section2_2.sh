#comments are for Rivanna (in progress)
wget https://raw.githubusercontent.com/guertinlab/Tn5bias/master/bedToOneEntryBed.py
#module load gcc/7.1.0
#module load seqoutbias/1.2.0
#module load bedtools/2.26.0
#module load intel/20.0  
#module load intelmpi/20.0
#module load python/3.7.7
#ijob -A guertinlab -c 1 -p standard -t 10:00:00

for i in C1_*_PE1_plus.bam
do
    name=$(echo $i | awk -F"_PE1_plus.bam" '{print $1}')
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
    seqOutBias hg38.fa ${name}_PE1_plus.bam --no-scale --out=no_scale.tbl --shift-counts \
                                 --bed=${name}_PE1_plus_shift_counts.bed \
                                 --bw=${name}_PE1_plus_shift_counts.bigWig --read-size=30
    seqOutBias hg38.fa ${name}_PE1_minus.bam --no-scale --out=no_scale.tbl --shift-counts \
                                 --bed=${name}_PE1_minus_shift_counts.bed \
                                 --bw=${name}_PE1_minus_shift_counts.bigWig --read-size=30
    seqOutBias hg38.fa ${name}_PE2_plus.bam --no-scale --out=no_scale.tbl --shift-counts \
                                 --bed=${name}_PE2_plus_shift_counts.bed \
                                 --bw=${name}_PE2_plus_shift_counts.bigWig --read-size=30
    seqOutBias hg38.fa ${name}_PE2_minus.bam --no-scale --out=no_scale.tbl --shift-counts \
                                 --bed=${name}_PE2_minus_shift_counts.bed \
                                 --bw=${name}_PE2_minus_shift_counts.bigWig --read-size=30
    python bedToOneEntryBed.py -i ${name}_PE1_plus_shift_counts.bed
    python bedToOneEntryBed.py -i ${name}_PE1_minus_shift_counts.bed
    python bedToOneEntryBed.py -i ${name}_PE2_plus_shift_counts.bed
    python bedToOneEntryBed.py -i ${name}_PE2_minus_shift_counts.bed
#these bed files can be used to get the sequence flanking ALL Tn5 insertion sites, 
#so we can see the precise nature of the sequence bias for each PE/strand combination
    awk '{$2 = $2 - 9; print}' ${name}_PE1_plus_shift_counts.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 20; print}' | grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_PE1_plus_shift_counts.fasta
    awk '{$2 = $2 - 9; print}' ${name}_PE1_minus_shift_counts.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 20; print}' |  grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_PE1_minus_not_scaled_shift_counts.fasta
    awk '{$2 = $2 - 9; print}' ${name}_PE2_plus_shift_counts.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 20; print}' | grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_PE2_plus_shift_counts.fasta
    awk '{$2 = $2 - 9; print}' ${name}_PE2_minus_shift_counts.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 20; print}' |  grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_PE2_minus_shift_counts.fasta
done


for i in C1_*_PE1_plus.bam
do
    name=$(echo $i | awk -F"_PE1_plus.bam" '{print $1}')
    echo $name
    seqOutBias hg38.fa ${name}_PE1_plus.bam --no-scale --out=no_scale.tbl \
                                 --bed=${name}_PE1_plus_no_shift.bed \
                                 --bw=${name}_PE1_plus_no_shift.bigWig --read-size=30
    seqOutBias hg38.fa ${name}_PE1_minus.bam --no-scale --out=no_scale.tbl \
                                 --bed=${name}_PE1_minus_no_shift.bed \
                                 --bw=${name}_PE1_minus_no_shift.bigWig --read-size=30
    seqOutBias hg38.fa ${name}_PE2_plus.bam --no-scale --out=no_scale.tbl \
                                 --bed=${name}_PE2_plus_no_shift.bed \
                                 --bw=${name}_PE2_plus_no_shift.bigWig --read-size=30
    seqOutBias hg38.fa ${name}_PE2_minus.bam --no-scale --out=no_scale.tbl \
                                 --bed=${name}_PE2_minus_no_shift.bed \
                                 --bw=${name}_PE2_minus_no_shift.bigWig --read-size=30
    python bedToOneEntryBed.py -i ${name}_PE1_plus_no_shift.bed
    python bedToOneEntryBed.py -i ${name}_PE1_minus_no_shift.bed
    python bedToOneEntryBed.py -i ${name}_PE2_plus_no_shift.bed
    python bedToOneEntryBed.py -i ${name}_PE2_minus_no_shift.bed
#these bed files can be used to get the sequence flanking ALL Tn5 insertion sites, 
#so we can see the precise nature of the sequence bias for each PE/strand combination
    awk '{$2 = $2 - 9; print}' ${name}_PE1_plus_no_shift.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 20; print}' | grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_PE1_plus_no_shift.fasta
    awk '{$2 = $2 - 9; print}' ${name}_PE1_minus_no_shift.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 20; print}' |  grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_PE1_minus_no_shift.fasta
    awk '{$2 = $2 - 9; print}' ${name}_PE2_plus_no_shift.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 20; print}' | grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_PE2_plus_no_shift.fasta
    awk '{$2 = $2 - 9; print}' ${name}_PE2_minus_no_shift.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 20; print}' |  grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_PE2_minus_no_shift.fasta
done


