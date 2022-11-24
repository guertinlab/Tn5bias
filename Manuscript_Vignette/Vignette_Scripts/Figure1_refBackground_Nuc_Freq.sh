#Create chrom.sizes file for hg38
faidx hg38.fa -i chromsizes > hg38.fa.chrom.sizes

#Find number of C in hg38:
grep -o "c\|C" hg38.fa | wc -l
#623727797
grep -o "c\|C" hg38.fa.chrom.sizes | wc -l
#455
#623727797-455 = 623727342

#Find number of G in hg38:
grep -o "g\|G" hg38.fa | wc -l
#626335226
grep -o "g\|G" hg38.fa.chrom.sizes | wc -l
#89
#626335226-89 = 626335137

#Find number of A in hg38:
grep -o "a\|A" hg38.fa | wc -l
#898285722
grep -o "a\|A" hg38.fa.chrom.sizes | wc -l
#303
#898285722-303 = 898285419

#Find number of T in hg38:
grep -o "t\|T" hg38.fa | wc -l
#900968146
grep -o "t\|T" hg38.fa.chrom.sizes | wc -l
#261
#900968146-261 = 900967885

#hg38
#623727342+626335137+898285419+900967885=3049315783
#C = 623727342/3049315783 = 0.2045467
#G = 626335137/3049315783 = 0.2054019
#A = 898285419/3049315783 = 0.2945859
#T = 900967885/3049315783 = 0.2954656
#0.2045467 + 0.2054019 + 0.2945859 + 0.2954656 = 1

#Create chrom.sizes file for hg38
faidx mm39.fa -i chromsizes > mm39.fa.chrom.sizes


#Find number of C in mm39:
grep -o "c\|C" mm39.fa | wc -l
#553008726
grep -o "c\|C" mm39.fa.chrom.sizes | wc -l
#61
#553008726-61 = 553008665

#Find number of G in mm39:
grep -o "g\|G" mm39.fa | wc -l
#553055984
grep -o "g\|G" mm39.fa.chrom.sizes | wc -l
#27
#553055984-27 = 553055957

#Find number of A in mm39:
grep -o "a\|A" mm39.fa | wc -l
#773810667
grep -o "a\|A" mm39.fa.chrom.sizes | wc -l
#18
#773810667-18 = 773810649

#Find number of T in mm39:
grep -o "t\|T" mm39.fa | wc -l
#774746512
grep -o "t\|T" mm39.fa.chrom.sizes | wc -l
#0
#774746512-0 = 774746512

#mm39
#553008665 + 553055957 + 773810649 + 774746512 = 2654621783
#C = 553008665/2654621783 = 0.2083192
#G = 553055957/2654621783 = 0.208337
#A = 773810649/2654621783 = 0.2914956
#T = 774746512/2654621783 = 0.2918482
#0.2083192 + 0.208337 + 0.2914956 + 0.2918482 = 1