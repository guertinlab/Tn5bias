
#if the next section works out, then it will eventually be incorporated into the loop above.

name=CEM-C7_untreated_rep2
#separate PE1 and PE2 reads
samtools view -b -f 0x0040 ${name}.bam > ${name}_PE1.bam
samtools view -b -f 0x0080 ${name}.bam > ${name}_PE2.bam


#separate plus and minus strand aligning reads:
samtools view -bh -F 20 ${name}_PE1.bam > ${name}_PE1_plus.bam
samtools view -bh -f 0x10 ${name}_PE1.bam > ${name}_PE1_minus.bam

samtools view -bh -F 20 ${name}_PE2.bam > ${name}_PE2_plus.bam
samtools view -bh -f 0x10 ${name}_PE2.bam > ${name}_PE2_minus.bam


#run seqOutBias to get the bed
seqOutBias hg38.fa ${name}_PE1_plus.bam --no-scale --bw=${name}_PE1_plus.bigWig --shift-counts --read-size=30
seqOutBias hg38.fa ${name}_PE1_minus.bam --no-scale --bw=${name}_PE1_minus.bigWig --shift-counts --read-size=30
seqOutBias hg38.fa ${name}_PE2_plus.bam --no-scale --bw=${name}_PE2_plus.bigWig --shift-counts --read-size=30
seqOutBias hg38.fa ${name}_PE2_minus.bam --no-scale --bw=${name}_PE2_minus.bigWig --shift-counts --read-size=30

# a python script to process the bed into a single entry per instance

wget https://raw.githubusercontent.com/guertinlab/Tn5bias/master/bedToOneEntryBed.py

python bedToOneEntryBed.py CEM-C7_untreated_rep2_PE1_plus_not_scaled.bed
python bedToOneEntryBed.py CEM-C7_untreated_rep2_PE1_minus_not_scaled.bed
python bedToOneEntryBed.py CEM-C7_untreated_rep2_PE2_plus_not_scaled.bed
python bedToOneEntryBed.py CEM-C7_untreated_rep2_PE2_minus_not_scaled.bed

#these bed files can be used to get the sequence flanking ALL Tn5 insertion sites, 
#so we can see the precise nature of the sequence bias for each PE/strand combination


tn5File=CEM-C7_untreated_rep2_PE1_plus_not_scaled.bed
awk '{$2 = $2 - 11; print}' $tn5File | awk '{OFS="\t";} {$3 = $2 + 20; print}' | grep -v - | \
    fastaFromBed -fi hg38.fa -s -bed stdin -fo CEM-C7_untreated_rep2_PE1_plus_not_scaled.fasta

tn5File=CEM-C7_untreated_rep2_PE1_minus_not_scaled.bed
awk '{$2 = $2 - 11; print}' $tn5File | awk '{OFS="\t";} {$3 = $2 + 20; print}' |  grep -v - | \
    fastaFromBed -fi hg38.fa -s -bed stdin -fo CEM-C7_untreated_rep2_PE1_minus_not_scaled.fasta


tn5File=CEM-C7_untreated_rep2_PE2_plus_not_scaled.bed
awk '{$2 = $2 - 11; print}' $tn5File | awk '{OFS="\t";} {$3 = $2 + 20; print}' | grep -v - | \
    fastaFromBed -fi hg38.fa -s -bed stdin -fo CEM-C7_untreated_rep2_PE2_plus_not_scaled.fasta

tn5File=CEM-C7_untreated_rep2_PE2_minus_not_scaled.bed
awk '{$2 = $2 - 11; print}' $tn5File | awk '{OFS="\t";} {$3 = $2 + 20; print}' |  grep -v - | \
    fastaFromBed -fi hg38.fa -s -bed stdin -fo CEM-C7_untreated_rep2_PE2_minus_not_scaled.fasta
