#This will make a .tbl file of all 3-mers in the genome
seqOutBias seqtable ../Figure1/hg38.fa --kmer-size=3 --plus-offset=3 --minus-offset=3 --read-size=76 --out=hg38.3.3.3.tbl
#This will dump all indexed 3-mer locations to a text file. This file will be LARGE
seqOutBias dump hg38.3.3.3.tbl > hg38.3.3.3.dump.txt
#This will take the plus and minus index for CAG and output a bed file for each
python dump_to_kmer_bed.py -i hg38.3.3.3.dump.txt -p 19 -m 31
#Combine the plus and minus bed files
cat hg38.3.3.3plus.19.bed hg38.3.3.3minus.31.bed > CAG_locations.bed
#Subset 400,000 random locations from combined bed files
shuf -n 400000 CAG_locations.bed > CAG_locations_rand400k.bed

#Run seqOutBias with masks designed to correct each peak (mask is shifted 4 base pairs):
seqOutBias ../Figure1/hg38.fa ../Figure1/C1_gDNA_rep1_plus.bam --custom-shift=4,-4 --strand-specific --kmer-mask=NNCN \
                                         --bed=C1_gDNA_rep1_plus_NNNXXXC.bed \
                                         --bw=C1_gDNA_rep1_plus_NNNXXXC.bigWig --read-size=76
seqOutBias ../Figure1/hg38.fa ../Figure1/C1_gDNA_rep1_minus.bam --custom-shift=4,-4 --strand-specific --kmer-mask=NNCN \
                                          --bed=C1_gDNA_rep1_minus_NNNXXXC.bed \
                                          --bw=C1_gDNA_rep1_minus_NNNXXXC.bigWig --read-size=76

seqOutBias ../Figure1/hg38.fa ../Figure1/C1_gDNA_rep1_plus.bam --custom-shift=4,-4 --strand-specific --kmer-mask=CXXXNNN \
                                         --bed=C1_gDNA_rep1_plus_NCNN.bed \
                                         --bw=C1_gDNA_rep1_plus_NCNN.bigWig --read-size=76
seqOutBias ../Figure1/hg38.fa ../Figure1/C1_gDNA_rep1_minus.bam --custom-shift=4,-4 --strand-specific --kmer-mask=CXXXNNN \
                                          --bed=C1_gDNA_rep1_minus_NCNN.bed \
                                          --bw=C1_gDNA_rep1_minus_NCNN.bigWig --read-size=76

seqOutBias ../Figure1/hg38.fa ../Figure1/C1_gDNA_rep1_plus.bam --custom-shift=4,-4 --strand-specific --kmer-mask=CXXXXXXXXNNN \
                                         --bed=C1_gDNA_rep1_plus_CXXXXNNN.bed \
                                         --bw=C1_gDNA_rep1_plus_CXXXXNNN.bigWig --read-size=76
seqOutBias ../Figure1/hg38.fa ../Figure1/C1_gDNA_rep1_minus.bam --custom-shift=4,-4 --strand-specific --kmer-mask=CXXXXXXXXNNN \
                                          --bed=C1_gDNA_rep1_minus_CXXXXNNN.bed \
                                          --bw=C1_gDNA_rep1_minus_CXXXXNNN.bigWig --read-size=76
