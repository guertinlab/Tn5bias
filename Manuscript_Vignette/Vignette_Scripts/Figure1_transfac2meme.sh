#Add 'XX' to end of each transfac document
endline='XX'
for i in *.transfac
do
    echo $endline >> $i
done

#Convert transfacs to memes
for i in *.transfac
do
    name=$(echo $i | awk -F".transfac" '{print $1}')
    echo $name
    transfac2meme $i > $name.meme

done