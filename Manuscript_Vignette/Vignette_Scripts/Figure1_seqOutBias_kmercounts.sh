setwd('EnzymeBias_masks')
#Determine all Tn5 scale factors
for i in Tn5*.tbl
do
    name=$(echo $i | awk -F".tbl" '{print $1}')
    echo $name
    seqOutBias table $i ../C1_gDNA_rep1.bam > ${name}_scale_factors.txt
done

#Determine all DNase scale factors
for i in DNase*.tbl
do
    name=$(echo $i | awk -F".tbl" '{print $1}')
    echo $name
    seqOutBias table $i ../DNase_Naked.bam > ${name}_scale_factors.txt
done

#Determine all Benzonase scale factors
for i in Benzonase*.tbl
do
    name=$(echo $i | awk -F".tbl" '{print $1}')
    echo $name
    seqOutBias table $i ../mm39_liver_Benzonase.bam > ${name}_scale_factors.txt
done

#Determine all Cyanase scale factors
for i in Cyanase*.tbl
do
    name=$(echo $i | awk -F".tbl" '{print $1}')
    echo $name
    seqOutBias table $i ../mm39_liver_Cyanase.bam > ${name}_scale_factors.txt
done

#Determine all MNase scale factors
for i in MNase*.tbl
do
    name=$(echo $i | awk -F".tbl" '{print $1}')
    echo $name
    seqOutBias table $i ../mm39_liver_MNase.bam > ${name}_scale_factors.txt
done