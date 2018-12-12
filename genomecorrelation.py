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
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr

# Adding local path to import program modules
sys.path.append(abspath(join(curdir, 'segmentanalysis')))
from segmentanalysis import segmentutils

chromosomeDataStr = 'Chromosome: {}, correlation: {:5.4f}, P-value: {:5.4f}'


def readEctopics(ectopicsFileName):
    ectopicsFile = open(ectopicsFileName)
    ectopics = {}
    for line in ectopicsFile:
        ID, value = line.strip().split()
        ectopics[ID] = float(value)
    return ectopics


def collectByChromosomes(intervalList):
    intervalsByChromosome = defaultdict(list)
    for interval in intervalList:
        intervalsByChromosome[interval.chromosome].append(interval)
    return intervalsByChromosome


# ANALYZING INPUT PARAMETERS
parser = argparse.ArgumentParser(description=programDescription, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("homologyFileName", type=str,
                    help="BED-formatted file holding homology data")
parser.add_argument("ectopicsFileName", type=str, nargs='?',
                    help="Tab-delimeted file containing ectopic contacts data")
args = parser.parse_args()

# Reading input data
homology = collectByChromosomes(segmentutils.readBedFile(args.homologyFileName))
ectopics = readEctopics(args.ectopicsFileName) if args.ectopicsFileName else None

# Converting to numpy arrays
nChromosomes = len(homology.keys())
fig, subplots = plt.subplots(nChromosomes, 1, squeeze=False)
subplots.shape = (nChromosomes,)

for i, (chromosome, intervalList) in enumerate(homology.items()):
    # Converting homolody metricks to NumPy arrays and calculating correlation
    homologyArr = np.array([float(interval.value) for interval in intervalList])
    # Plotting data
    currPlot = subplots[i]
    title = 'Chromosome: ' + chromosome

    positions = np.array([(interval.start + interval.stop) / 2 for interval in intervalList])
    # Plotting homology data
    color = 'tab:blue'
    currPlot.set_xlabel('Position (b.p.)')
    currPlot.set_ylabel('Homology', color=color)
    currPlot.plot(positions, homologyArr, color=color)
    currPlot.tick_params(axis='y', labelcolor=color)
    if ectopics is not None:  # Plotting ectopics
        ectopicsArr = np.array([ectopics[interval.ID] for interval in intervalList if interval.ID in ectopics.keys()])
        if len(ectopicsArr) > 0:
            corrCoef, pValue = pearsonr(homologyArr, ectopicsArr)
            print(chromosomeDataStr.format(chromosome, corrCoef, pValue))
            title = chromosomeDataStr.format(chromosome, corrCoef, pValue)
            ax2 = currPlot.twinx()  # instantiate a second axes that shares the same x-axis
            color = 'tab:red'
            ax2.set_ylabel('ectopics', color=color)  # we already handled the x-label with ax1
            ax2.plot(positions, ectopicsArr, color=color)
            ax2.tick_params(axis='y', labelcolor=color)
            fig.tight_layout()  # otherwise the right y-label is slightly clipped
    currPlot.set_title(title)
plt.show()
