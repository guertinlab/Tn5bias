for ba in C1*
do
	name=$(echo ${ba:0:-4})
	#Set window size for FASTA	
	winsize=21
	minusshift=10	
	echo $name
	echo Window size is $winsize\bp
	echo Window begins $minusshift\bp before insertion event
	
	#separate PE1 and PE2 reads
	samtools view -b -f 0x0040 ${name}.bam > ${name}_PE1.bam
	samtools view -b -f 0x0080 ${name}.bam > ${name}_PE2.bam


	#separate plus and minus strand aligning reads:
	samtools view -bh -F 20 ${name}_PE1.bam > ${name}_PE1_plus.bam
	samtools view -bh -f 0x10 ${name}_PE1.bam > ${name}_PE1_minus.bam

	samtools view -bh -F 20 ${name}_PE2.bam > ${name}_PE2_plus.bam
	samtools view -bh -f 0x10 ${name}_PE2.bam > ${name}_PE2_minus.bam


	#run seqOutBias to get the bed
	seqOutBias hg38.fa ${name}_PE1_plus.bam --no-scale --bw=${name}_PE1_plus.bigWig --read-size=30
	seqOutBias hg38.fa ${name}_PE1_minus.bam --no-scale --bw=${name}_PE1_minus.bigWig --read-size=30
	seqOutBias hg38.fa ${name}_PE2_plus.bam --no-scale --bw=${name}_PE2_plus.bigWig --read-size=30
	seqOutBias hg38.fa ${name}_PE2_minus.bam --no-scale --bw=${name}_PE2_minus.bigWig --read-size=30



	python bedToOneEntryBed.py -i ${name}_PE1_plus_not_scaled.bed
	python bedToOneEntryBed.py -i ${name}_PE1_minus_not_scaled.bed
	python bedToOneEntryBed.py -i ${name}_PE2_plus_not_scaled.bed
	python bedToOneEntryBed.py -i ${name}_PE2_minus_not_scaled.bed

	#these bed files can be used to get the sequence flanking ALL Tn5 insertion sites, 
	#so we can see the precise nature of the sequence bias for each PE/strand combination


	tn5File=${name}_PE1_plus_not_scaled.oneentry.bed
	awk -v m=$minusshift '{$2 = $2 - m; print}' $tn5File | \
		awk -v w=$winsize '{OFS="\t";} {$3 = $2 + w; print}' | grep -v - | \
    		fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_PE1_plus_not_scaled.fasta

	tn5File=${name}_PE1_minus_not_scaled.oneentry.bed
	awk -v m=$minusshift '{$2 = $2 - m; print}' $tn5File | \
		awk -v w=$winsize '{OFS="\t";} {$3 = $2 + w; print}' | grep -v - | \
    		fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_PE1_minus_not_scaled.fasta


	tn5File=${name}_PE2_plus_not_scaled.oneentry.bed
	awk -v m=$minusshift '{$2 = $2 - m; print}' $tn5File | \
		awk -v w=$winsize '{OFS="\t";} {$3 = $2 + w; print}' | grep -v - | \
    		fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_PE2_plus_not_scaled.fasta

	tn5File=${name}_PE2_minus_not_scaled.oneentry.bed
	awk -v m=$minusshift '{$2 = $2 - m; print}' $tn5File | \
		awk -v w=$winsize '{OFS="\t";} {$3 = $2 + w; print}' | grep -v - | \
    		fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_PE2_minus_not_scaled.fasta
	
	#reverse complement the negative strands for comparison with positive strands
    	
	awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${name}_PE1_minus_not_scaled.fasta | \
        	fastx_reverse_complement -o ${name}_PE1_minus_no_shift_RC.fasta
    	
	awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${name}_PE2_minus_not_scaled.fasta | \
        	fastx_reverse_complement -o ${name}_PE2_minus_no_shift_RC.fasta

done
