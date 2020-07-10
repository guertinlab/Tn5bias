for ba in C1*
do
	name=$(echo $ba)
	#Set window size for FASTA	
	winsize=20
	minusshift=$(expr $winsize / 2 - 1)	
	echo $name
	echo ${name:0:-4}
	echo Window size is $winsize\bp
	echo Window begins $minusshift\bp before insertion event
	
	#separate PE1 and PE2 reads
	samtools view -b -f 0x0040 ${name:0:-4}.bam > ${name:0:-4}_PE1.bam
	samtools view -b -f 0x0080 ${name:0:-4}.bam > ${name:0:-4}_PE2.bam


	#separate plus and minus strand aligning reads:
	samtools view -bh -F 20 ${name:0:-4}_PE1.bam > ${name:0:-4}_PE1_plus.bam
	samtools view -bh -f 0x10 ${name:0:-4}_PE1.bam > ${name:0:-4}_PE1_minus.bam

	samtools view -bh -F 20 ${name:0:-4}_PE2.bam > ${name:0:-4}_PE2_plus.bam
	samtools view -bh -f 0x10 ${name:0:-4}_PE2.bam > ${name:0:-4}_PE2_minus.bam


	#run seqOutBias to get the bed
	seqOutBias hg38.fa ${name:0:-4}_PE1_plus.bam --no-scale --bw=${name:0:-4}_PE1_plus.bigWig --shift-counts --read-size=30
	seqOutBias hg38.fa ${name:0:-4}_PE1_minus.bam --no-scale --bw=${name:0:-4}_PE1_minus.bigWig --shift-counts --read-size=30
	seqOutBias hg38.fa ${name:0:-4}_PE2_plus.bam --no-scale --bw=${name:0:-4}_PE2_plus.bigWig --shift-counts --read-size=30
	seqOutBias hg38.fa ${name:0:-4}_PE2_minus.bam --no-scale --bw=${name:0:-4}_PE2_minus.bigWig --shift-counts --read-size=30



	python bedToOneEntryBed.py -i ${name:0:-4}_PE1_plus_not_scaled.bed
	python bedToOneEntryBed.py -i ${name:0:-4}_PE1_minus_not_scaled.bed
	python bedToOneEntryBed.py -i ${name:0:-4}_PE2_plus_not_scaled.bed
	python bedToOneEntryBed.py -i ${name:0:-4}_PE2_minus_not_scaled.bed

	#these bed files can be used to get the sequence flanking ALL Tn5 insertion sites, 
	#so we can see the precise nature of the sequence bias for each PE/strand combination


	tn5File=${name:0:-4}_PE1_plus_not_scaled.oneentry.bed
	awk -v m=$minusshift '{$2 = $2 - m; print}' $tn5File | awk -v w=$winsize '{OFS="\t";} {$3 = $2 + w; print}' | grep -v - | \
    		fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name:0:-4}_PE1_plus_not_scaled.fasta

	tn5File=${name:0:-4}_PE1_minus_not_scaled.oneentry.bed
	awk -v m=$minusshift '{$2 = $2 - m; print}' $tn5File | awk -v w=$winsize '{OFS="\t";} {$3 = $2 + w; print}' | grep -v - | \
    		fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name:0:-4}_PE1_minus_not_scaled.fasta


	tn5File=${name:0:-4}_PE2_plus_not_scaled.oneentry.bed
	awk -v m=$minusshift '{$2 = $2 - m; print}' $tn5File | awk -v w=$winsize '{OFS="\t";} {$3 = $2 + w; print}' | grep -v - | \
    		fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name:0:-4}_PE2_plus_not_scaled.fasta

	tn5File=${name:0:-4}_PE2_minus_not_scaled.oneentry.bed
	awk -v m=$minusshift '{$2 = $2 - m; print}' $tn5File | awk -v w=$winsize '{OFS="\t";} {$3 = $2 + w; print}' | grep -v - | \
    		fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name:0:-4}_PE2_minus_not_scaled.fasta



done
