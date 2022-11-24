fimo --thresh 0.00005 --text Tn5_sepcat_bias.meme \
                             hg38.fa > Tn5_sepcat_bias_transfac_00005_fimo.txt

fimo --thresh 0.0001 --text DNase_sepcat_bias.meme \
                            hg38.fa > DNase_sepcat_bias_transfac_0001_fimo.txt

fimo --thresh 0.00007 --text Cyanase_sepcat_bias.meme \
                             mm39.fa > Cyanase_sepcat_bias_transfac_00007_fimo.txt

fimo --thresh 0.00005 --text Benzonase_sepcat_bias.meme \
                             mm39.fa > Benzonase_sepcat_bias_transfac_00005_fimo.txt

fimo --thresh 0.0001 --text MNase_sepcat_bias.meme \
                            mm39.fa > MNase_sepcat_bias_transfac_0001_fimo.txt
