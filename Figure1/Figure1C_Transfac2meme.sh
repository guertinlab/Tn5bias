endline='XX'
for i in *.transfac
do
    echo $endline >> $i

done

for i in *.transfac
do
    name=$(echo $i | awk -F".transfac" '{print $1}')
    echo $name
    dm transfac2meme $i > $name.meme

done