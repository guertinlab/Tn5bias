#Make a fasta file of sequences 50bp upstream of all plus strand cutsites
awk '{$2 = $2 - 50; print}' ../Figure1/C1_gDNA_rep1_plus_not_scaled.oneentry.bed | \
awk '{OFS="\t";} {$3 = $2 + 3; print}' | grep -v - | \
fastaFromBed -fi ../Figure1/hg38.fa -s -bed stdin -fo C1_gDNA_rep1_50bp_plus.fasta

#Make a fasta file of sequences 50bp upstream of all minus strand cutsites
awk '{$2 = $2 - 50; print}' ../Figure1/C1_gDNA_rep1_minus_not_scaled.oneentry.bed | \
awk '{OFS="\t";} {$3 = $2 + 3; print}' |  grep -v - | \
fastaFromBed -fi ../Figure1/hg38.fa -s -bed stdin -fo C1_gDNA_rep1_50bp_minus.fasta

#Reverse complement all minus strand sequences
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' C1_gDNA_rep1_50bp_minus.fasta | \
fastx_reverse_complement -o C1_gDNA_rep1_50bp_minus_RC.fasta

#Concatenate the plus and minus strand fasta files together
cat C1_gDNA_rep1_50bp_minus_RC.fasta C1_gDNA_rep1_50bp_plus.fasta > C1_gDNA_rep1_50bp_sepcat.fasta



################################################################################################################
#Make a fasta file of sequences 100bp upstream of all plus strand cutsites
awk '{$2 = $2 - 100; print}' ../Figure1/C1_gDNA_rep1_plus_not_scaled.oneentry.bed | \
awk '{OFS="\t";} {$3 = $2 + 3; print}' | grep -v - | \
fastaFromBed -fi ../Figure1/hg38.fa -s -bed stdin -fo C1_gDNA_rep1_100bp_plus.fasta

#Make a fasta file of sequences 100bp upstream of all minus strand cutsites
awk '{$2 = $2 - 100; print}' ../Figure1/C1_gDNA_rep1_minus_not_scaled.oneentry.bed | \
awk '{OFS="\t";} {$3 = $2 + 3; print}' |  grep -v - | \
fastaFromBed -fi ../Figure1/hg38.fa -s -bed stdin -fo C1_gDNA_rep1_100bp_minus.fasta

#Reverse complement all minus strand sequences
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' C1_gDNA_rep1_100bp_minus.fasta | \
fastx_reverse_complement -o C1_gDNA_rep1_100bp_minus_RC.fasta

#Concatenate the plus and minus strand fasta files together
cat C1_gDNA_rep1_100bp_minus_RC.fasta C1_gDNA_rep1_100bp_plus.fasta > C1_gDNA_rep1_100bp_sepcat.fasta