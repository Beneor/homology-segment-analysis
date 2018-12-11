#!/usr/bin/env python3

programDescription = '''
This script produces a set of short (5 - 30 by default) random fragments for the specified chromosome segment,
splits chromosome into a set of long (10 kb by default) chunks,
seeks random segments within chromosome parts,
computes frequencies of occurence for segments of different length within chromosome parts,
makes table for the normalized frequencies,
and calculates Pearson correlation between fragment frequencies for chromosome discs and ectopic contact frequencies
'''

import argparse
import sys
from os.path import join, abspath, curdir

import numpy as np
import matplotlib.pyplot as plt 
from scipy.stats.stats import pearsonr

# Adding local path to import program modules
sys.path.append(abspath(join(curdir, 'segmentanalysis')))
from segmentanalysis import segmentutils


def readEctopics(ectopicsFileName):
    ectopicsFile = open(ectopicsFileName)
    ectopics = {}
    for line in ectopicsFile:
        ID, value = line.strip().split()
        ectopics[ID] = float(value)
    return ectopics


# 0. ANALYZING INPUT PARAMETERS
parser = argparse.ArgumentParser(description=programDescription, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("homologyFileName", type=str,
                    help="FASTA file (may be gzipped) containing genome sequence")
parser.add_argument("ectopicsFileName", type=str,
                    help="segment of chromosome to analyze (X:16113516:16900779) in BED-file notation (starting from 0, end is not included)")

args = parser.parse_args()

ectopics = readEctopics(args.ectopicsFileName)
homology = segmentutils.readBedFile(args.homologyFileName)

# Converting to numpy arrays

homologyArr = np.array([float(interval.value) for interval in homology])
ectopicsArr = np.array([ectopics[interval.ID]*1000 for interval in homology])
coor = np.array([(interval.start+interval.stop)/2 for interval in homology])

corrCoef, pValue = pearsonr(homologyArr, ectopicsArr )

plt.plot(coor, homologyArr, coor, ectopicsArr )
plt.show()


print("Correlation coefficient: ", corrCoef)
print("P-Value: ", pValue)
