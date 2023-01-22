weblogo -f DNase_sepcat_bias.transfac -D transfac -o DNase_sepcat_bias.eps \
        --color-scheme classic -F eps -S 0.35 -Y YES -s large \
        --composition "{'A':29.5, 'C':20.5, 'G':20.5, 'T':29.5}" --logo-font Arial-BoldMT \
        --annotate '-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15' \
        --number-fontsize 14 --rotate-numbers YES --fineprint '' --ylabel '' --errorbars NO

weblogo -f Tn5_sepcat_bias.transfac -D transfac -o Tn5_sepcat_bias.eps \
        --color-scheme classic -F eps -S 0.35 -Y YES -s large \
        --composition "{'A':29.5, 'C':20.5, 'G':20.5, 'T':29.5}" --logo-font Arial-BoldMT \
        --annotate '-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15' \
        --number-fontsize 14 --rotate-numbers YES --fineprint '' --ylabel '' --errorbars NO

weblogo -f Benzonase_sepcat_bias.transfac -D transfac -o Benzonase_sepcat_bias.eps \
        --color-scheme classic -F eps -S 0.35 -Y YES -s large \
        --composition "{'A':29.2, 'C':20.8, 'G':20.8, 'T':29.2}" --logo-font Arial-BoldMT \
        --annotate '-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15' \
        --number-fontsize 14 --rotate-numbers YES --fineprint '' --ylabel '' --errorbars NO

weblogo -f Cyanase_sepcat_bias.transfac -D transfac -o Cyanase_sepcat_bias.eps \
        --color-scheme classic -F eps -S 0.35 -Y YES -s large \
        --composition "{'A':29.2, 'C':20.8, 'G':20.8, 'T':29.2}" --logo-font Arial-BoldMT \
        --annotate '-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15' \
        --number-fontsize 14 --rotate-numbers YES --fineprint '' --ylabel '' --errorbars NO

weblogo -f MNase_sepcat_bias.transfac -D transfac -o MNase_sepcat_bias.eps \
        --color-scheme classic -F eps -S 0.35 -Y YES -s large \
        --composition "{'A':29.2, 'C':20.8, 'G':20.8, 'T':29.2}" --logo-font Arial-BoldMT \
        --annotate '-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15' \
        --number-fontsize 14 --rotate-numbers YES --fineprint '' --ylabel '' --errorbars NO
################################################################################################################
weblogo -f Tn5_sepcat_bias.transfac -D transfac -o Tn5_sepcat_bias_greatestval.eps \
        --color-scheme classic -F eps -S 0.31 -Y YES -s large \
        --composition "{'A':29.5, 'C':20.5, 'G':20.5, 'T':29.5}" --logo-font Arial-BoldMT \
        --annotate '-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15' \
        --number-fontsize 14 --rotate-numbers YES --fineprint '' --ylabel '' --errorbars NO

weblogo -f Benzonase_sepcat_bias.transfac -D transfac \
        -o Benzonase_sepcat_bias_greatestval.eps \
        --color-scheme classic -F eps -S 0.16 -Y YES -s large \
        --composition "{'A':29.2, 'C':20.8, 'G':20.8, 'T':29.2}" --logo-font Arial-BoldMT \
        --annotate '-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15' \
        --number-fontsize 14 --rotate-numbers YES --fineprint '' --ylabel '' --errorbars NO

weblogo -f Cyanase_sepcat_bias.transfac -D transfac \
        -o Cyanase_sepcat_bias_greatestval.eps --color-scheme classic \
        -F eps -S 0.13 -Y YES -s large \
        --composition "{'A':29.2, 'C':20.8, 'G':20.8, 'T':29.2}" --logo-font Arial-BoldMT \
        --annotate '-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15' \
        --number-fontsize 14 --rotate-numbers YES --fineprint '' --ylabel '' --errorbars NO

weblogo -f MNase_sepcat_bias.transfac -D transfac \
        -o MNase_sepcat_bias_greatestval.eps --color-scheme classic \
        -F eps -S 0.31 -Y YES -s large  \
        --composition "{'A':29.2, 'C':20.8, 'G':20.8, 'T':29.2}" --logo-font Arial-BoldMT \
        --annotate '-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15' \
        --number-fontsize 14 --rotate-numbers YES --fineprint '' --ylabel '' --errorbars NO

################################################################################################################
weblogo -f DNase_sepcat_bias.transfac -D transfac \
        -o DNase_sepcat_bias_EQUIPROBABLE_BACKGROUND.eps --color-scheme classic \
        -F eps -S 0.31 -Y YES -s large --composition equiprobable --logo-font Arial-BoldMT \
        --annotate '-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15' \
        --number-fontsize 14 --rotate-numbers YES --fineprint '' --ylabel '' --errorbars NO

weblogo -f Tn5_sepcat_bias.transfac -D transfac \
        -o Tn5_sepcat_bias_EQUIPROBABLE_BACKGROUND.eps --color-scheme classic \
        -F eps -S 0.18 -Y YES -s large --composition equiprobable --logo-font Arial-BoldMT \
        --annotate '-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15' \
        --number-fontsize 14 --rotate-numbers YES --fineprint '' --ylabel '' --errorbars NO

weblogo -f Benzonase_sepcat_bias.transfac -D transfac \
        -o Benzonase_sepcat_bias_EQUIPROBABLE_BACKGROUND.eps --color-scheme classic \
        -F eps -S 0.1 -Y YES -s large --composition equiprobable --logo-font Arial-BoldMT \
        --annotate '-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15' \
        --number-fontsize 14 --rotate-numbers YES --fineprint '' --ylabel '' --errorbars NO

weblogo -f Cyanase_sepcat_bias.transfac -D transfac \
        -o Cyanase_sepcat_bias_EQUIPROBABLE_BACKGROUND.eps --color-scheme classic \
        -F eps -S 0.1 -Y YES -s large --composition equiprobable --logo-font Arial-BoldMT \
        --annotate '-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15' \
        --number-fontsize 14 --rotate-numbers YES --fineprint '' --ylabel '' --errorbars NO

weblogo -f MNase_sepcat_bias.transfac -D transfac \
        -o MNase_sepcat_bias_EQUIPROBABLE_BACKGROUND.eps --color-scheme classic \
        -F eps -S 0.42 -Y YES -s large --composition equiprobable --logo-font Arial-BoldMT \
        --annotate '-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15' \
        --number-fontsize 14 --rotate-numbers YES --fineprint '' --ylabel '' --errorbars NO
