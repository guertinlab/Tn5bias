#Make a new directory for rule ensemble masks
mkdir RE_masks
while IFS= read -r line; do
    echo $line
		cp 'Tn5_'$line'.bw' RE_masks
done < RulesEnsemble_top10_masks.txt
cd RE_masks
#convert bigWigs into bedGraphs
for bw in Tn5*.bw
do
	name=$(echo $bw | awk -F".bw" '{print $1}')
  echo $name
  bigWigToBedGraph $name'.bw' $name'.bedGraph'
done
#Make a text file with all bedgraph names and order
for bw in *.bedGraph
do
	name=$(echo $bw | awk -F".bedGraph" '{print $1}')
	echo $name
	echo $name >> Tn5_top10_names.txt
done
#Combine all bedGraphs into a single file
bedtools unionbedg -i *.bedGraph > Tn5_top10_union.bedGraph