#!/usr/bin/env python3

import analysismodules.analysisfunctions3 as an
import os, shutil
import glob
import copy

#Parameters  #!!! Change these parameters, if necessary
frd='../examples/Cytocountsdata'
ectd='../examples/Xectopicsdata'
lengthstuple1 = (10,15,20,25,30,35,40,45,50,55,60) # Fragment lengths (full sample)
lengthstuple2 = (25,30) # Test sample
strain1 = 'CS'
strain2 = 'agn'
d='discs' #X chromosome divided into 119 discs (excluding disc 20B, == disc 202)
z='zones' #X chromosome divided into 19 zonew (excluding zone 20)
al='all' #All fragments
nr='nr' #No repeats
s='specific'
u='unspecific'
pardicts=[{"0":strain1,"1":strain2},{"0":d,"1":z},{"0":al,"1":nr},{"0":s,"1":u}]  #Dictionary containing alternative parameters for statistical analysis 


# 2.3. Mann-Whitney U analysis of correlations in the same strains for different fragment lengths for nearest discs
resultdir='./MW-analysis/nearest'
os.mkdir(resultdir) 
for t in an.parcombmdiff(n=4,dic=pardicts,v=1):
    if 'zones' in t:
        continue    #Excluding data for zones from analysis
    if 'nr' in t:
        continue
    for n in lengthstuple2:
        an.nearestcompare(*t,n,x='Nearest')
    for f in glob.glob('*'+'.xls'):
        shutil.move(f,resultdir)
    
