#! /sw/bin/python
import re
import string
import sys
import getopt
import os
import itertools
import glob

def function(fqFileName):
    infile=open(fqFileName, 'r')
    outfile=open(str.split(fqFileName,'.bed')[0]+'.oneentry.bed', 'w')
    while 1:
        line=infile.readline()
        if not line: break
        for i in range(int(line.split()[4])):
            outfile.write('%s'%(line))
    infile.close()
    outfile.close()
    
    
def main(argv):
    try:
        opts, args = getopt.getopt(argv, "i:h", ["infile=","help"])
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)
    infile = False
    for opt, arg in opts:
        if opt in ('-i', '--infile'):
            infile = arg
        elif opt in ('-h', '--help'):
            print('\n./bedToOneEntryBed.py -i test_PE1_plus_not_scaled.bed')
            sys.exit()
    if infile:
        print(infile)
        function(infile)
if __name__ == "__main__":
    main(sys.argv[1:])
