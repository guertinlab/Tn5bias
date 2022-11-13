#! /sw/bin/python
import re
import string
import sys
import getopt
import os
import itertools
import glob

def function(dumpFileName, kmerp, kmerm):
    infile=open(dumpFileName, 'r')
    outfileplus=open(str.split(dumpFileName,'.dump')[0]+ 'plus.' + kmerp + '.bed', 'w')
    outfileminus=open(str.split(dumpFileName,'.dump')[0]+ 'minus.' + kmerm + '.bed', 'w')
    sze = int(infile.readline().split()[2])
    plusoff = infile.readline().split()[2]
    minusoff = int(infile.readline().split()[2])
    readlen = int(infile.readline().split()[2])
    while 1:
        line=infile.readline()
        if not line: break
        splitline = line.split()
        if line.startswith('>'):
            chr = splitline[0].split('>')[1]
        else:
            if splitline[1] == kmerp:
                start = int(splitline[0]) - int(plusoff)
                outfileplus.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(chr, str(start), str(start + sze), splitline[1], str(kmerp), '+'))
#still workin gon this. do I need to RC on my own?
            if splitline[2] == kmerm:
                start = int(splitline[0]) + readlen
                outfileminus.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(chr, str(start), str(start + sze), splitline[2], str(kmerm), '-'))
    infile.close()
    outfileplus.close()
    outfileminus.close()
    
def main(argv):
    try:
        opts, args = getopt.getopt(argv, "i:p:m:h", ["dump=", "kmerplus=","kmerminus=","help"])
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)
    infile = False
    for opt, arg in opts:
        if opt in ('-i', '--dump'):
            infile = arg
        elif opt in ('-p', '--kmerplus'):
            inf2 = arg
        elif opt in ('-m', '--kmerminus'):
            inf3 = arg
        elif opt in ('-h', '--help'):
            print('\n python dump_to_kmer_bed.py -i hg38.3.3.3.dump.test.txt -kp 64 -km 1')
            sys.exit()
    if infile and inf2 and inf3:
        print(infile)
        function(infile, inf2, inf3)
if __name__ == "__main__":
    main(sys.argv[1:])
