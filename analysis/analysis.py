#!/usr/bin/env python3

import analysismodules.analysisfunctions as an
import os, shutil
import glob
import copy

programDescription = '''
! Place this script into the folder 'analysis' (if absent) and run it in that folder using Python3. 

This script performes analysis of correlations between frequencies of matching fragments (FMF) and frequencies of ectopic contacts (FEC).
'''

#Parameters  #!!! Change these parameters, if necessary
# Aho-corasick is run by shell script
frd='../examples/Cytocountsdata'  # Frequency of mathced fragments across genome - need to fill in before starting analysis. Populated by Filesmodification.py - 

ectd='../examples/Xectopicsdata'  # Input folder with ectopics
lengthstuple1 = (10,15,20,25,30,35,40,45,50,55,60) # Fragment lengths (full sample)
lengthstuple2 = (55,60) # Test sample
strain1 = 'CS'    #names of strains for Xectopicsdata
strain2 = 'agn'   #



d='discs' #X chromosome divided into 119 discs (excluding disc 20B, == disc 202)  # In cytocounts and XEctopics
z='zones' #X chromosome divided into 19 zonew (excluding zone 20)  # In cytocounts and XEctopics

al='all' #All fragments # folders for fragments total 
nr='nr' # and No repeats

s='specific'    # Type of correlations to list - for the same band 
u='unspecific'  # For not same bandes
pardicts=[{"0":strain1,"1":strain2},{"0":d,"1":z},{"0":al,"1":nr},{"0":s,"1":u}]  #Dictionary containing alternative parameters for statistical analysis 


# 1. Spearman correlation analysis 
#1.1. Correlations for fragments frequencies and ectopic contacts frequencies (all X chromosome)
os.mkdir('./Correlations/')
for t in an.parcomb(n=4,dic=pardicts):
    an.correlation(*t,lengths=lengthstuple1,fragmentsdir=frd,ectopicsdir=ectd)
    
# 1.2. Mann-Whitney statistical analysis of difference between correlation values and probabilities

os.mkdir('./MW-analysis/')
# 1.2.1. Mann-Whitney U analysis of correlations for the data sets with the same chief parameters (different fragment lengths) 
resultdir='./MW-analysis/diff-lengths'
os.mkdir(resultdir)
for t in an.parcomb(n=4,dic=pardicts):
    if 'zones' in t: #Excluding data for zones from analysis
        continue
    an.statmw(*t,*t,x='Correlations',u=True)
for f in glob.glob('*'+'.xls'):
    an.shutil.move(f,resultdir)
    
# 1.2.2. Mann-Whitney U analysis of correlations for the same fragment lengths with one different chief parameter.
resultdir='./MW-analysis/diff-parameters'
os.mkdir(resultdir)
for t in an.parcombmdiff(n=4,dic=pardicts,v=1):
    if 'zones' in t:
        continue    #Excluding data for zones from analysis
    an.statmw(*t,x='Correlations',u=False)
for f in glob.glob('*'+'.xls'):
    shutil.move(f,resultdir)
    
    
# 2.1. Correlations for fragments frequencies and ectopic contacts frequencies (nearby discs within the different areas)
os.mkdir('./Nearest/')

for t in an.parcomb(n=4,dic=pardicts):
    if 'zones' in t:
        continue
    if 'specific' in t:
        an.neardisccor (*t,lengths=lengthstuple1,fragmentsdir=frd,ectopicsdir=ectd)
    if 'unspecific' in t:
        an.neardiscunspecifcor(*t,lengths=lengthstuple1,fragmentsdir=frd,ectopicsdir=ectd)
    
# 2.2. Mann-Whitney U analysis of correlations in the same strains for different fragment lengths for nearest discs
resultdir='./MW-analysis/nearest'
os.mkdir(resultdir) 
for t in an.parcombmdiff(n=4,dic=pardicts,v=1):
    if 'zones' in t:
        continue    #Excluding data for zones from analysis
    for n in lengthstuple1:
        an.nearestcompare(*t,n,x='Nearest')
    for f in glob.glob('*'+'.xls'):
        shutil.move(f,resultdir)
        
