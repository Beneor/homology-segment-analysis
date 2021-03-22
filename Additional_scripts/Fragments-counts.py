#!/usr/bin/env python3

import glob
import os

programDescription = '''
! Place this script into the folder containing the fragments you want to analyze (e.g. 'l50-fragments') and run it in that folder using Python3. 

This script produces a list of fragments making impact into FMF.
For each unique fragment, the script output its sequence and the number of occurence in descending order ( >10). The result file is in csv format and can be opened by LibreOffice Calc tables 
'''

filesuffix = '.txt' 
filemask = '*'+ filesuffix
d = ","
files = glob.glob(filemask)

r = open(r'Result.csv','w')
a={}
b={}
c={}
l=[]
for filename in files:
    file1 = open(filename,'r')
    for line in file1:
        x = str(line).strip()
        if x == 0:
            continue
        y = x.count(d)+1 #The number of fragments separated by comma for a given position
        x1=x[x.find('\t')+1:x.rfind('\t')] #Fragment sequence
        if x1 in a:
            y+=a[x1]
        a[x1]=y #The number of fragment sequences in dictionary with key x1 (fragment sequence)
    file1.close()
for i in a:
    if a[i]>10:
        b[i] = a[i]

for i in b:
    j=b[i]
    c[j]=i
    l.append(j)
s=set(l)
l=list(s)
l.sort(reverse=True)

for i in l:
    r.write(str(i))
    r.write('\t')
    r.write(str(c[i]))
    r.write('\n')
r.close()
