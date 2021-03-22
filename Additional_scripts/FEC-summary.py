#!/usr/bin/env python3
import glob
import os, shutil

filesuffix= '.csv'
filemask='*'+filesuffix
files=glob.glob(filemask)
resultfile = open (r'Average-ECF.csv2','w')

programDescription = '''
! Place this script into the folder containing files with frequencies of ectopic contacts (FEC) you want to analyze (e.g. 'Xectopicsdata/CS/discs') and run it in that folder using Python3.

This script generates the summary FEC for each band.
'''

c={}
e=[]
for filename in files:
    f=open(filename, "r")
    fn = int(filename.strip(filesuffix))
    a = []
    for line in f:
        x = str(line)
        ecf=float(x[x.find('\t'):])
        a.append(ecf)
    c[fn]=a
    e.append(fn)
    f.close()
e.sort()
b=[]
for i in e:
   b.append(sum(c[i]))
for i in b:
    resultfile.write (str(i))
    resultfile.write('\n')
resultfile.close()

