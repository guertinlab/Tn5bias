####FIMO transcription factor motifs using hg38
for meme in ../Figure4/*.meme
do
    name=$(echo $meme | awk -F".meme" '{print $1}' | awk -F"../Figure4/" '{print $2}')
    echo $name
    fimo --thresh 0.001 --text ${meme} ../Figure1/hg38.fa > ${name}_FIMO.txt
done
#Remove DUX4
rm DUX4_FIMO.txt
##wget repeatmasker coordinates
wget  https://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz
gunzip -c hg38.fa.out.gz > hg38_repeatmasker.fa
awk  '{ OFS="\t" } {print $5, $6, $7, $11}' hg38_repeatmasker.fa | tail +4 > hg38_repeatmasker.bed
rm hg38_repeatmasker.fa
###Remove repeatmasker regions from all FIMO files
for file in *_FIMO.txt
do
  TF=$(echo $file | awk -F"_FIMO.txt" '{print $1}')
  echo ${TF}
  awk  '{ OFS="\t" } {print $3, $4, $5, $2, $1, $6, $7, $8, $9}' $file | tail +2 > ${TF}_FIMO.bed
  bedtools intersect -wa -v \
  -a ${TF}_FIMO.bed \
  -b hg38_repeatmasker.bed \
   > ${TF}.removed.bed
  echo ${TF}.removed.bed
  awk  '{ OFS="\t" } {print $5, $4, $1, $2, $3, $6, $7, $8, "", $9}' ${TF}.removed.bed > ${TF}_nohead_FIMO.bed
  head -n 1 AR_fimo.txt | cat - ${TF}_nohead_FIMO.bed > ${TF}_rm_FIMO.txt
  rm ${TF}_FIMO.bed ${TF}.removed.bed ${TF}_nohead_FIMO.bed
done
