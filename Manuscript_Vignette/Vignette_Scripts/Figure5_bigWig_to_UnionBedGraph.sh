#convert bigWigs into bedGraphs
for bw in ../Figure1/DNase*.bigWig
do
	name=$(echo $bw | awk -F".bigWig" '{print $1}' | awk -F"../Figure1/" '{print $2}')
  echo $name
  bigWigToBedGraph '../Figure1/'$name'.bigWig' $name'.bedGraph'
done

#Make a text file with all bedgraph names and order
for bw in *.bedGraph
do
	name=$(echo $bw | awk -F".bedGraph" '{print $1}')
	echo $name
	echo $name >> DNase_allmasks_names.txt
done

#Combine all bedGraphs into a single file
bedtools unionbedg -i *.bedGraph > DNase_5mer_union.bedGraph