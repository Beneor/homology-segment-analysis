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
from scipy.stats.stats import pearsonr, spearmanr

# Adding local path to import program modules
sys.path.append(abspath(join(curdir, 'segmentanalysis')))
from segmentanalysis import segmentutils

chromosomeDataStr = 'Chromosome: {}, {} correlation: {:5.4f}, P-value: {:5.4f}'

correlationMethods = {'Spearman': spearmanr,
                      'Pearson': pearsonr
                      }


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

parser.add_argument("-r", "--correlation", type=str, choices=correlationMethods.keys(), default='Spearman',
                    help="Method for correlation calculations")
parser.add_argument("-c", "--chromosomes", type=str,
                    help="List of chromosomes divided by , to view (view all if not set. Example:2L,X")
parser.add_argument("-l", "--labels", action='store_true',
                    help="Put labels from homology data on X axis")
parser.add_argument("-m", "--mark", action='store_true',
                    help="Add marks to graph")
parser.add_argument("-x", "--xinterval", type=int, default=1000,
                    help="Interval to put ticks on X-axis, in kilobases")
parser.add_argument("-s", "--savefig", type=str,
                    help="file name to save resulting picture")

args = parser.parse_args()

xinterval = args.xinterval * 1000

# Reading input data
homology = collectByChromosomes(segmentutils.readBedFile(args.homologyFileName))
# Filtering by provided set of chromosomes,

if args.chromosomes is not None:
    chrSet = set(args.chromosomes.split(','))
    homology = {key: value for key, value in homology.items() if key in chrSet}

ectopics = readEctopics(args.ectopicsFileName) if args.ectopicsFileName else None

# Grating subgraphs
nChromosomes = len(homology.keys())
fig, subplots = plt.subplots(nChromosomes, 1, squeeze=False)
subplots.shape = (nChromosomes,)  # Resizing to linear array

for i, (chromosome, intervalList) in enumerate(homology.items()):
    # Plotting data
    currPlot = subplots[i]
    title = 'Chromosome: ' + chromosome
    color = 'tab:blue'
    xLabel = 'Position (Mbases)'
    currPlot.set_ylabel('Homology', color=color)

    # Converting to NumPy and plotting homology data
    positions = np.array([(interval.start + interval.stop) / 2 for interval in intervalList])
    homologyArr = np.array([float(interval.value) for interval in intervalList])
    currPlot.plot(positions, homologyArr, marker='o' if args.mark else None, color=color, label='Homology level')

    # Setting graph labels, parameters et.c
    plt.sca(currPlot)  # Selecting current plot
    plt.tick_params(axis='y', labelcolor=color)
    if args.labels:  # Setting labels from interval IDs
        labels = [interval.ID for interval in intervalList]
        plt.xticks(positions, labels, fontsize='xx-small', rotation='vertical')
        xLabel = 'Intervals'
    else:  # Setting labels to megabases
        locs = np.arange(0, intervalList[-1].stop, xinterval)
        plt.xticks(locs, map(lambda x: "{:2.1f}".format(x), locs / 1e6))
    currPlot.set_xlabel(xLabel)
    plt.xlim([0, intervalList[-1].stop])

    if ectopics is not None:  # Plotting ectopics
        ectopicsArr = np.array([ectopics[interval.ID] for interval in intervalList if interval.ID in ectopics.keys()])
        if len(ectopicsArr) > 0:
            corrCoef, pValue = correlationMethods[args.correlation](homologyArr, ectopicsArr)
            title = chromosomeDataStr.format(chromosome, args.correlation, corrCoef, pValue)
            print(title)
            ax2 = currPlot.twinx()  # instantiate a second axes that shares the same x-axis
            color = 'tab:red'
            ax2.set_ylabel('ectopics', color=color)  # we already handled the x-label with ax1
            ax2.plot(positions, ectopicsArr, marker='x' if args.mark else None, color=color,
                     label='Ectopic contacts frequency')
            ax2.tick_params(axis='y', labelcolor=color)
            fig.tight_layout()  # otherwise the right y-label is slightly clipped
    currPlot.set_title(title)
    #currPlot.legend()
fig.tight_layout()
if args.savefig:
    plt.savefig(args.savefig)
plt.show()
