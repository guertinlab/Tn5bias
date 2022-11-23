#Download Benzonase, Cyanase, MNase mm39 (mouse) data
fasterq-dump SRR535737
fasterq-dump SRR535738
fasterq-dump SRR535739
fasterq-dump SRR535740
fasterq-dump SRR535741
fasterq-dump SRR535742
fasterq-dump SRR535743
fasterq-dump SRR535744
fasterq-dump SRR5723785

#Rename Benzonase, Cyanase, MNase data from SRR number
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
cat *Benz* > mm39_liver_Benzonase.fastq
cat *Cyan* > mm39_liver_Cyanase.fastq
cat *MNase* > mm39_liver_MNase.fastq

#Zip each file
gzip *ase.fastq

#build index genome for mm39
bowtie2-build mm39.fa mm39

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
    seqOutBias mm39.fa ${name}.bam --no-scale --shift-counts \
                                --strand-specific --bed=${name}.bed \
                                 --bw=${name}.bigWig --read-size=76
    seqOutBias mm39.fa ${name}_plus.bam --no-scale --shift-counts \
                                --strand-specific --bed=${name}_plus.bed \
                                 --bw=${name}_plus.bigWig --read-size=76
    seqOutBias mm39.fa ${name}_minus.bam --no-scale --shift-counts \
                                --strand-specific --bed=${name}_minus.bed \
                                 --bw=${name}_minus.bigWig --read-size=76
#Make each read its own entry in a bed file
#e.g. 2 of the same read becomes 2 entries of the same read
    python bedToOneEntryBed.py -i ${name}_not_scaled.bed
    python bedToOneEntryBed.py -i ${name}_plus_not_scaled.bed
    python bedToOneEntryBed.py -i ${name}_minus_not_scaled.bed
#These bed files can be used to get the sequence flanking ALL Tn5 insertion sites,
#so we can see the precise nature of the sequence bias for each PE/strand combination
#we shift the window coordinatess for each bed file to accomodate the 31bp for the figure
    awk '{$2 = $2 - 10; print}' ${name}_not_scaled.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 21; print}' | grep -v - | \
       fastaFromBed -fi mm39.fa -s -bed stdin -fo ${name}.fasta
    awk '{$2 = $2 - 10; print}' ${name}_plus_not_scaled.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 21; print}' | grep -v - | \
       fastaFromBed -fi mm39.fa -s -bed stdin -fo ${name}_plus.fasta
    awk '{$2 = $2 - 11; print}' ${name}_minus_not_scaled.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 21; print}' |  grep -v - | \
       fastaFromBed -fi mm39.fa -s -bed stdin -fo ${name}_minus.fasta
#First convert all sequences in the minus fasta into uppercase letters,
#then reverse complement all minus strand sequences.
    awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${name}_minus.fasta | \
       fastx_reverse_complement -o ${name}_minus_RC.fasta
#Concatenate the plus and minus strand fasta files together
    cat ${name}_minus_RC.fasta ${name}_plus.fasta > ${name}_sepcat.fasta
done

#Download Tn5 dataset
fasterq-dump SRR5123141

#Zip Tn5 data
gzip *fastq

#Rename Tn5 data
mv SRR5123141_1.fastq.gz C1_gDNA_rep1_PE1.fastq.gz
mv SRR5123141_2.fastq.gz C1_gDNA_rep1_PE2.fastq.gz


#build index genome
bowtie2-build hg38.fa hg38
bowtie2-build chrM.fa chrM

#First align Tn5 data to mitochondrial chromosome, then sort by chr/start
#Convert bam file to fastq file, then align to hg38 genome
for fq in *PE1.fastq.gz
do
	name=$(echo $fq | awk -F"_PE1.fastq.gz" '{print $1}')
	echo $name
	bowtie2 -p 6 -x chrM -1 ${name}_PE1.fastq.gz -2 ${name}_PE2.fastq.gz \
		| samtools view -b - | samtools sort - -o $name.sorted.chrM.bam
	samtools index $name.sorted.chrM.bam
	samtools view -b $name.sorted.chrM.bam '*' | samtools sort -n - \
		| bamToFastq -i - -fq ${name}_PE1.chrM.fastq -fq2 ${name}_PE2.chrM.fastq
	gzip *.chrM.fastq
	rm $name.sorted.chrM.bam
	bowtie2 -p 6 --maxins 500 -x hg38 -1 ${name}_PE1.chrM.fastq.gz -2 ${name}_PE2.chrM.fastq.gz \
		| samtools view -bS - | samtools sort -n - | samtools fixmate -m - - \
		| samtools sort - | samtools markdup -r - ${name}.bam

done

#First, split data into plus and minus aligned reads
#Run seqOutBias with no scale and shift of 4,-4 to get each read centered
#on cut site. Then process same as other enzymes above.
for i in C1*.bam
do
    name=$(echo $i | awk -F".bam" '{print $1}')
    echo $name


		samtools view -bh -F 20 ${name}.bam > ${name}_plus.bam
    samtools view -bh -f 0x10 ${name}.bam > ${name}_minus.bam
    seqOutBias hg38.fa ${name}.bam --no-scale --custom-shift=4,-4 \
                                --strand-specific --bed=${name}.bed \
                                 --bw=${name}.bigWig --read-size=76
    seqOutBias hg38.fa ${name}_plus.bam --no-scale --custom-shift=4,-4 \
                                --strand-specific --bed=${name}_plus.bed \
                                 --bw=${name}_plus.bigWig --read-size=76
    seqOutBias hg38.fa ${name}_minus.bam --no-scale --custom-shift=4,-4 \
                                --strand-specific --bed=${name}_minus.bed \
                                 --bw=${name}_minus.bigWig --read-size=76
    python bedToOneEntryBed.py -i ${name}.bed
		python bedToOneEntryBed.py -i ${name}_plus.bed
    python bedToOneEntryBed.py -i ${name}_minus.bed
#we shift the window coordinatess for each bed file to accomodate the 31bp for the figure
    awk '{$2 = $2 - 15; print}' ${name}_not_scaled.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 31; print}' | grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}.fasta
		awk '{$2 = $2 - 15; print}' ${name}_plus_not_scaled.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 31; print}' | grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_plus.fasta
    awk '{$2 = $2 - 15; print}' ${name}_minus_not_scaled.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 31; print}' |  grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_minus.fasta

		awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${name}_minus.fasta | \
	     fastx_reverse_complement -o ${name}_minus_RC.fasta
#Concatenate the plus and minus strand fasta files together
	  cat ${name}_RC.fasta ${name}_plus.fasta > ${name}_sepcat.fasta

done

#Download naked DNase data
fasterq-dump SRR769954

#Rename DNase data from SRR number
mv SRR769954.fastq DNase_Naked.fastq

#Zip
gzip DNase_Naked.fastq

#Align each dataset to hg38 genome
for fq in DNase_*.fastq.gz
do
    name=$(echo $fq | awk -F".fastq.gz" '{print $1}')
    echo $name
    bowtie2 -p 6 --maxins 500 -x hg38 -U ${name}.fastq.gz  \
       | samtools view -bS - | samtools sort - -o $name.bam
done

#First, split each dataset into plus and minus aligned reads.
#Then, run seqOutBias on plus/minus and unseparated strands to get
##the bedfile which has chromosome coordinates for each read.
#From the bedfile, make each read a unique entry using bedToOneEntryBed.py
#Using awk commands, shift each position back 10bp then forward 21 in order to
##make a 'window' around the cut site.
#In order to get each sequence, we then convert from fasta to bed format for each data set
for i in DNase_*.bam
do
    name=$(echo $i | awk -F".bam" '{print $1}')
    echo $name
    samtools view -bh -F 20 ${name}.bam > ${name}_plus.bam
    samtools view -bh -f 0x10 ${name}.bam > ${name}_minus.bam
    seqOutBias hg38.fa ${name}.bam --no-scale --shift-counts \
                                --strand-specific --bed=${name}_unscaled.bed \
                                 --bw=${name}_unscaled.bigWig --read-size=76
    seqOutBias hg38.fa ${name}_plus.bam --no-scale --shift-counts \
                                --strand-specific --bed=${name}_plus_unscaled.bed \
                                 --bw=${name}_plus_unscaled.bigWig --read-size=76
    seqOutBias hg38.fa ${name}_minus.bam --no-scale --shift-counts \
                                --strand-specific --bed=${name}_minus_unscaled.bed \
                                 --bw=${name}_minus_unscaled.bigWig --read-size=76
#Make each read its own entry in a bed file
#e.g. 2 of the same read becomes 2 entries of the same read
    python bedToOneEntryBed.py -i ${name}_not_scaled.bed
    python bedToOneEntryBed.py -i ${name}_plus_not_scaled.bed
    python bedToOneEntryBed.py -i ${name}_minus_not_scaled.bed
#these bed files can be used to get the sequence flanking ALL DNase insertion sites,
#so we can see the precise nature of the sequence bias for each PE/strand combination
#we shift the window coordinatess for each bed file to accomodate the 31bp for the figure
    awk '{$2 = $2 - 15; print}' ${name}_not_scaled.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 31; print}' | grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}.fasta
    awk '{$2 = $2 - 15; print}' ${name}_plus_not_scaled.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 31; print}' | grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_plus.fasta
    awk '{$2 = $2 - 16; print}' ${name}_minus_not_scaled.oneentry.bed | \
       awk '{OFS="\t";} {$3 = $2 + 31; print}' |  grep -v - | \
       fastaFromBed -fi hg38.fa -s -bed stdin -fo ${name}_minus.fasta
#First convert all sequences in the minus fasta into uppercase letters,
#then reverse complement all minus strand sequences.
    awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${name}_minus.fasta | \
       fastx_reverse_complement -o ${name}_minus_RC.fasta
#Concatenate the plus and minus strand fasta files together
    cat ${name}_minus_RC.fasta ${name}_plus.fasta > ${name}_sepcat.fasta
done
