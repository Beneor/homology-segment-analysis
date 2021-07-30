#!/usr/bin/env python3

import segmentanalysis.analysisfunctions as an
import os, shutil
import glob
import argparse

programDescription = '''
This script performes analysis of correlations between frequencies of matching fragments (FMF) and frequencies of ectopic contacts (FEC).
It require the pre-calculated frequencies of matching fragments (FMF), which can be calculated by the Filesmodification.py script. 
The Aho-corasick runs for Filesmodification.py should be prepared by the shell script
'''

# 0. ANALYZING INPUT PARAMETERS
parser = argparse.ArgumentParser(description=programDescription, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("fmfGenome", type=str,
                    help="Frequency of matching fragments (FMF)", default="./examples/Cytocountsdata", nargs='?')
parser.add_argument("ectopicsData", type=str,
                    help="Input folder with ectopics data", default="./examples/Xectopicsdata", nargs='?')

#parser.add_argument("outDir", type=str, help="Directory to write results WARNING: Currently not implemented. Please don't use", default=".", nargs='?')
parser.add_argument("strain1", type=str, help="First strain to analyse. Must be the same as in folder with ectopics data", default="CS", nargs='?')
parser.add_argument("strain2", type=str, help="Second strain to analyse. Must be the same as in folder with ectopics data", default="agn", nargs='?')

parser.add_argument("-s", "--fragmentsizes", type=str, default='30,50',
                    help='Set of fragment sizes to search')

parser.add_argument("-N",  "--nearest", action='store_true',
                    help='Perform nearest calculations. WARNING: in the current version is very time-consuming due to non-optimal algorithm.')

args = parser.parse_args()
args.outDir = '.' # Quick fix to workaround not implemented output directory selection
frd = args.fmfGenome # fragmentsdir
ectd = args.ectopicsData  # ectopicsdir

lengthstuple1 = tuple([int(size) for size in args.fragmentsizes.split(',')])

strain1 = args.strain1  # names of strains for Xectopicsdata
strain2 = args.strain2  #

d = 'discs'  # X chromosome divided into 119 discs (excluding disc 20B, == disc 202)  # In cytocounts and XEctopics
z = 'zones'  # X chromosome divided into 19 zonew (excluding zone 20)  # In cytocounts and XEctopics

al = 'all'  # All fragments # folders Но я могу for fragments total
nr = 'nr'  # and No repeats

s = 'specific'  # Type of correlations to list - for the same band
u = 'unspecific'  # For not same bandes

pardicts = [{"0": strain1, "1": strain2}, {"0": d, "1": z}, {"0": al, "1": nr},
            {"0": s, "1": u}]  # Dictionary containing alternative parameters for statistical analysis

# 1. Spearman correlation analysis
# 1.1. Correlations for fragments frequencies and ectopic contacts frequencies (all X chromosome)

os.mkdir(os.path.join(args.outDir, 'Correlations'))

for t in an.parcomb(n=4, dic=pardicts):
    an.correlation(*t, lengths=lengthstuple1, fragmentsdir=frd, ectopicsdir=ectd)

# 1.2. Mann-Whitney statistical analysis of difference between correlation values and probabilities

os.mkdir(os.path.join(args.outDir, 'MW-analysis'))
# 1.2.1. Mann-Whitney U analysis of correlations for the data sets with the same chief parameters (different fragment lengths) 
resultdir = os.path.join(args.outDir, 'MW-analysis', 'diff-lengths')
os.mkdir(resultdir)
print("Analysing MW")
for t in an.parcomb(n=4, dic=pardicts):
    #print("Analyzing combination: ",t)
    
    if not os.path.exists(os.path.join(frd,t[1],t[2])):
        print("Skipping analysis for non-existing combination " + '-'.join(t))
        continue
    print("Analyzing combination " + '-'.join(t))

    an.statmw(*t, *t, x='Correlations', u=True)
for f in glob.glob('*' + '.xls'):
    an.shutil.move(f, resultdir)

# 1.2.2. Mann-Whitney U analysis of correlations for the same fragment lengths with one different chief parameter.
resultdir = './MW-analysis/diff-parameters'
os.mkdir(resultdir)
print("Analysing MW- diff-parameters")
for t in an.parcombmdiff(n=4, dic=pardicts, v=1):
    
    if not os.path.exists(os.path.join(frd,t[1],t[2])) or not os.path.exists(os.path.join(frd,t[5],t[6])):
        print("Skipping analysis for non-existing combination " + '-'.join(t[:4]) + " against " + '-'.join(t[4:]))
        continue
    print("Analyzing combination "  + '-'.join(t[:4]) + " against " + '-'.join(t[4:]))

    an.statmw(*t, x='Correlations', u=False)

for f in glob.glob('*' + '.xls'):
    shutil.move(f, resultdir)

# 2.1. Correlations for fragments frequencies and ectopic contacts frequencies (nearby discs within the different areas)
if args.nearest:
    os.mkdir('./Nearest/')

    for t in an.parcomb(n=4, dic=pardicts):
        if not os.path.exists(os.path.join(frd,t[1],t[2])):
            print("Skipping 'nearest' analysis for non-existing combination " + '-'.join(t))
            continue
        print("Analyzing 'nearest' combination " + '-'.join(t))
    
        if 'specific' in t:
            an.neardisccor(*t, lengths=lengthstuple1, fragmentsdir=frd, ectopicsdir=ectd)
        if 'unspecific' in t:
            an.neardiscunspecifcor(*t, lengths=lengthstuple1, fragmentsdir=frd, ectopicsdir=ectd)

    # 2.2. Mann-Whitney U analysis of correlations in the same strains for different fragment lengths for nearest discs
    resultdir = './MW-analysis/nearest'
    os.mkdir(resultdir)
    for t in an.parcombmdiff(n=4, dic=pardicts, v=1):
        if 'zones' in t:
            continue  # Excluding data for zones from analysis
        for n in lengthstuple1:
            an.nearestcompare(*t, n, x='Nearest')
        for f in glob.glob('*' + '.xls'):
            shutil.move(f, resultdir)
