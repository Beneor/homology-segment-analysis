#!/usr/bin/env python3

import argparse
import sys
import glob
import re

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr, spearmanr

from os.path import join, abspath, curdir

# Adding local path to import program modules
sys.path.append(abspath(join(curdir, "segmentanalysis")))
from segmentanalysis import utils, fragmentstatistics

programDescription = """

See details in program description in README.md
"""


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


def findSubstringWithOverlaps(sub, string):
    pos = 0
    indexes = []
    while True:
        pos = string.find(sub, pos)
        if pos > -1:
            indexes.append(pos)
            pos += 1
        else:
            break
    return indexes
    # return [m.start() for m in re.finditer(sub, string, overlapped=True)]


# 0. ANALYZING INPUT PARAMETERS
parser = argparse.ArgumentParser(
    description=programDescription, formatter_class=argparse.RawDescriptionHelpFormatter
)
# parser.add_argument("dataFolder", type=str,
# help="Output data forlder for segmentanalysis run")
parser.add_argument(
    "fastaFileName",
    type=str,
    help="FASTA file (may be gzipped) containing genome sequence",
)
parser.add_argument(
    "fragment", type=str, help="sequence of fragment to find correlations"
)
parser.add_argument("cytomap", type=str, help="BED file containing cytobands")
parser.add_argument(
    "-v",
    "--verbose",
    action="store_true",
    help="Print additional information to stdout",
)
parser.add_argument(
    "-c",
    "--chunk",
    type=float,
    default=10.0,
    help="Chunk size to divide chromosome, in kilobases",
)
parser.add_argument(
    "-e",
    "--exclude",
    type=str,
    help="Coordinated of interval where we need to exclude signal, start:stop",
)

args = parser.parse_args()

# I. Loading data
if args.verbose:
    print("Loading genome from FASTA file")
genome = segmentutils.readFasta(segmentutils.openMaybeGzipped(args.fastaFileName))
chunkSize = int(args.chunk * 1000)


# II. Searching for fragment locations

chrPositionsChunks = {}
for chromosome in genome:
    # Now we search only for single fragment - no need to use Aho-Corasick search
    fragmentPositionsFor = np.fromiter(
        findSubstringWithOverlaps(args.fragment, genome[chromosome]), dtype=np.int64
    )
    fragmentPositionsRev = np.fromiter(
        findSubstringWithOverlaps(
            segmentutils.revcomp(args.fragment), genome[chromosome]
        ),
        dtype=np.int64,
    )
    if args.verbose:
        lFor = len(fragmentPositionsFor)
        lRev = len(fragmentPositionsRev)
        print(
            "Chromosome {} (forward) - {} occurences, ({}-forward, {}-reverse)".format(
                chromosome, lFor + lRev, lFor, lRev
            )
        )
    # print(fragmentPositionsFor)
    # open(args.fragment+'-'+chromosome+'-for.txt','w').write(np.array2string(fragmentPositionsFor))

    fragmentPositions = np.concatenate((fragmentPositionsFor, fragmentPositionsRev))
    fragmentPositions.sort()
    chrPositionsChunks[chromosome] = (
        fragmentPositions // chunkSize
    )  # Dividing each element by chunk size using numpy syntax

# III. Converting to counts
counts = segmentstatistics.chunksToCounts(chrPositionsChunks, genome, chunkSize)
if args.exclude is not None:
    exclude = segmentutils.strToBed(args.exclude, separator=":")
    segmentstatistics.excludeIntervalFromCounts(counts, exclude, chunkSize)
    if args.verbose:
        print(
            "Setting to zero counts for chunks {}:{}-{}".format(
                exclude.chromosome,
                exclude.start // chunkSize,
                exclude.stop // chunkSize,
            )
        )
        print(
            "Corrected total fragments number for chromosome {}: {}".format(
                exclude.chromosome, np.sum(counts[exclude.chromosome])
            )
        )

# III. Converting counts to cytobands
cytomap = segmentutils.readBedFile(args.cytomap)
cytoCounts = segmentstatistics.countsToCytomap(cytomap, counts, chunkSize)
if args.verbose:
    print("Total matches: ", np.sum(cytoCounts))
cytoCountsFileName = "{}.{}.cytocounts.txt".format(args.fastaFileName, args.fragment)
segmentutils.dumpCytoCouns(cytoCountsFileName, cytomap, cytoCounts)
