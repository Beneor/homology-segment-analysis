#!/usr/bin/env python3

import argparse
import sys
from os.path import join, abspath, curdir
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr, spearmanr

# Adding local path to import program modules
sys.path.append(abspath(join(curdir, 'segmentanalysis')))
from segmentanalysis import utils


programDescription = '''
This script makes plot of normalisez matching frequencies calculated by segmentanalysis.py script. 
Script is able to make plots for counts normalized by chunks or by cytobands. 
If some experiments data for genome (for example ectopic contacts frequency) are avaliable, 
script plots experimental data on the same grach and calculates correlation between fragments matching frequencies
adn experimental data.

See details in program description in README.md
'''

chromosomeDataStr = 'Chromosome: {}, {} correlation: {:5.4f}, P-value: {:5.4f}'
correlationMethods = {'Spearman': spearmanr,
                      'Pearson': pearsonr
                      }

def readEctopics(ectopicsFileName):
    """
    Reads ectopics or nay similar genome data from tab-delimited file
    :param ectopicsFileName: file name to read data
    :return: dictionary of regions named adn regions data
    """
    ectopicsFile = open(ectopicsFileName)
    ectopics = {}
    for line in ectopicsFile:
        ID, value = line.strip().split()
        ectopics[ID] = float(value)
    return ectopics


def collectByChromosomes(intervalList):
    """
    Grouping list of genome intervals by chromosomes
    :param intervalList: list of BED intervals, read from file
    :return: dictionary:
        key: chromosome name
        value: list of intervals for particular chromosome
    """
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

# Sixe of interval for x-axis ticks
xinterval = args.xinterval * 1000

# Reading input data
homology = collectByChromosomes(segmentutils.readBedFile(args.homologyFileName))
# Filtering by provided set of chromosomes,

if args.chromosomes is not None:
    chrSet = set(args.chromosomes.split(','))
    homology = {key: value for key, value in homology.items() if key in chrSet}

ectopics = readEctopics(args.ectopicsFileName) if args.ectopicsFileName else None
#print(args.homologyFileName)
#print(args.ectopicsFileName)

# Grating subgraphs
nChromosomes = len(homology.keys())
fig, subplots = plt.subplots(nChromosomes, 1, squeeze=False)
subplots.shape = (nChromosomes,)  # Resizing to linear array

for i, (chromosome, intervalList) in enumerate(homology.items()):
    # Plotting data
    currPlot = subplots[i]
    title = 'Chromosome: ' + chromosome # Default title for graphs
    color = 'tab:blue'
    currPlot.set_ylabel('Homology', color=color)

    # Converting to NumPy and plotting homology data
    positions = np.array([(interval.start + interval.stop) / 2 for interval in intervalList])
    homologyArr = np.array([float(interval.value) for interval in intervalList])
    currPlot.plot(positions, homologyArr, marker='o' if args.mark else None, color=color, label='Homology level')

    # Setting graph labels, parameters et c.
    plt.sca(currPlot)  # Selecting current plot
    plt.tick_params(axis='y', labelcolor=color)
    if args.labels:  # Setting labels from interval IDs
        labels = [interval.ID for interval in intervalList]
        plt.xticks(positions, labels, fontsize='xx-small', rotation='vertical')
        xLabel = 'Intervals'
    else:  # Setting labels to megabases
        locs = np.arange(0, intervalList[-1].stop, xinterval)
        plt.xticks(locs, map(lambda x: "{:2.1f}".format(x), locs / 1e6))
        xLabel = 'Position (Mbases)'
    currPlot.set_xlabel(xLabel)
    plt.xlim([0, intervalList[-1].stop])

    if ectopics is not None:  # Plotting ectopics
        ectopicsArr = np.array([ectopics[interval.ID] for interval in intervalList
                                if interval.ID in ectopics.keys()]) # Converting to numpy array
        if len(ectopicsArr) > 0:
            corrCoef, pValue = correlationMethods[args.correlation](homologyArr, ectopicsArr)
            title = chromosomeDataStr.format(chromosome, args.correlation, corrCoef, pValue) # Changing title
            print(title)
            ax2 = currPlot.twinx()  # instantiate a second axes that shares the same x-axis
            color = 'tab:red'
            ax2.set_ylabel('ectopics', color=color)  # we already handled the x-label with ax1
            ax2.plot(positions, ectopicsArr, marker='x' if args.mark else None, color=color,
                     label='Ectopic contacts frequency')
            ax2.tick_params(axis='y', labelcolor=color)
            fig.tight_layout()  # otherwise the right y-label is slightly clipped
    currPlot.set_title(title)
    # currPlot.legend() #TODO: find how to add global legend for all graphs
fig.tight_layout()
if args.savefig:
    plt.savefig(args.savefig)
else:
    plt.show()

