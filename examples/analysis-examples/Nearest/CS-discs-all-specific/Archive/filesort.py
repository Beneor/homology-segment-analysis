#!/usr/bin/env python3
import glob
import os, shutil

lengthstuple1 = (10,15,20,25,30,35,40,45,50,55,60)

for i in lengthstuple1:
    j=str(i)
    os.mkdir(j)
    filesuffix= j+'-nearest.csv2'
    filemask='*'+filesuffix
    files=glob.glob(filemask)
    for f in glob.glob('*'+filesuffix):
        shutil.move(f,j)
