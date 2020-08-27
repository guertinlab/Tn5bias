#build index genome
bowtie2-build hg38.fa hg38
bowtie2-build chrM.fa chrM

for fq in *PE1.fastq.gz
do
	name=$(echo $fq | awk -F"_PE1.fastq.gz" '{print $1}')
	echo $name
	fastp -f 4 -F 4 -i ${name}_PE1.fastq.gz -I ${name}_PE2.fastq.gz \
		-o ${name}_PE1.trimmed.fastq.gz -O ${name}_PE2.trimmed.fastq.gz
	bowtie2 -x chrM -1 ${name}_PE1.trimmed.fastq.gz -2 ${name}_PE2.trimmed.fastq.gz \
		| samtools view -b - | samtools sort - -o $name.sorted.chrM.bam
	samtools index $name.sorted.chrM.bam
	samtools view -b $name.sorted.chrM.bam '*' | samtools sort -n - \
		| bamToFastq -i - -fq ${name}_PE1.chrM.fastq -fq2 ${name}_PE2.chrM.fastq
	gzip *.chrM.fastq
	rm $name.sorted.chrM.bam
	bowtie2 --maxins 500 -x hg38 -1 ${name}_PE1.chrM.fastq.gz -2 ${name}_PE2.chrM.fastq.gz \
		| samtools view -bS - | samtools sort -n - | samtools fixmate -m - - \
		| samtools sort - | samtools markdup -r - $name.bam
	seqOutBias hg38.fa $name.bam --no-scale --bw=${name}.bigWig --shift-counts --read-size=30
done
