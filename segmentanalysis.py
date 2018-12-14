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
import gzip
import os
from os.path import join, abspath, curdir

# Adding local path to import program modules
sys.path.append(abspath(join(curdir, 'segmentanalysis')))
from segmentanalysis import segmentutils
from segmentanalysis import segmentsearch, segmentstatistics

# 0. ANALYZING INPUT PARAMETERS

parser = argparse.ArgumentParser(description=programDescription, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("fastaFileName", type=str,
                    help="FASTA file (may be gzipped) containing genome sequence")
parser.add_argument("segment", type=str,
                    help="segment of chromosome to analyze (X:16113516:16900779) in BED-file notation (starting from 0, end is not included)")
parser.add_argument("cytomap", type=str,nargs='?',
                    help="BED file containing cytobands")

parser.add_argument("-v", "--verbose", action='store_true',
                    help='Print additional information to stdout')
parser.add_argument("-s", "--fragmentsizes", type=str, default='5,10,15,20,25,30',
                    help='Set of fragment sizes to search')
parser.add_argument("-d", "--fragmentdensity", type=int, default=5,
                    help='Average frequency of fragments in letters. Default 5 means that in average each 5-th letter will be start of fragment')
parser.add_argument("-c", "--chunk", type=float, default=10.0,
                    help='Chunk size to divide chromosome, in kilobases')
parser.add_argument("-m", "--mergedirections", action='store_true', default=True,
                    help='Output combined ncounts by both directions')
parser.add_argument("--dump", action='store_true', default=True,
                    help='Save found fragments and positions to files (used only for fragments >= mindumpsize)')
parser.add_argument("--mindumpsize", type=int, default=10,
                    help='Minumum size of fragment to store in file')

args = parser.parse_args()

# I. INPUT FILES PREPARATION

# 0. Reading genome sequence from FASTA[.gz] file
if args.fastaFileName.endswith('.gz'):
    fastaFile = gzip.open(args.fastaFileName, 'rt')
else:
    fastaFile = open(args.fastaFileName, 'r')

if args.verbose:
    print('Loading genome from FASTA file')
genome = segmentutils.readFasta(fastaFile)
fastaFile.close()

# 1. Selection of chromosome segment
segment = segmentutils.segmentStrToGGenomeInterval(args.segment, genome)

segmentSeqFor = genome[segment.chromosome][segment.start:segment.stop]
segmentSeqs = {'dir': segmentSeqFor, 'rev': segmentutils.revcomp(segmentSeqFor)}

# 2. Size of chunks and selected fragments
fragmentSizes = [int(size) for size in args.fragmentsizes.split(',')]
chunkSize = int(args.chunk * 1000)

# 3. Folder to save results
outputFolder = '{}.{}-{}-{}'.format(args.fastaFileName, segment.chromosome, segment.start, segment.stop)
if os.path.exists(outputFolder):
    print('Output folder {} already exists\nPlase delete and resubmit'.format(outputFolder))
    exit(0)
os.makedirs(outputFolder)

if args.cytomap is not None:
    cytomap = segmentutils.readBedFile(args.cytomap)
    if args.verbose:
        print('Grouping counts by cytomap regions')


# II. ITERATION THROUGH ALL FRAGMENT SIZES AND DIRECTIONS

# Generating random fragments from selected region
for fragmentSize in fragmentSizes:
    normalizedCounts = {}
    cytomapCounts = {}
    for direction, segmentSeq in segmentSeqs.items():
        # III. GENERATING FRAGMENTS
        if args.verbose:
            print('Generating fragments in {} direction for size: {}'.format(direction, fragmentSize))
        fragments = segmentsearch.chooseFragments(segmentSeq, fragmentSize, args.fragmentdensity)

        # IV. SEARCH OF MATCHING FRAGMENTS
        chrFragmentsPositions = segmentsearch.searchFragments(genome, fragments, args.verbose)
        if args.dump and fragmentSize >= args.mindumpsize:
            if args.verbose:
                print('Dumping fragments to text file')
            fragmentsFileName = '{}/fragments.l{:02d}-{}.txt'.format(outputFolder, fragmentSize, direction)
            segmentutils.dumpFragmentsToFile(fragmentsFileName, chrFragmentsPositions)
        # v. CONVEERTING TO CHUNKS AND NORMALIZING
        chrPositionsChunks = segmentstatistics.locationsToChunks(chrFragmentsPositions, chunkSize)
        counts = segmentstatistics.chunksToCounts(chrPositionsChunks, genome, chunkSize)
        # Removing counts inside fragment
        if args.verbose:
            print('Setting to zero counts for chunks', segment.chromosome,
                  ':', segment.start // chunkSize, '-', segment.stop // chunkSize)
        for chunk in range(segment.start // chunkSize, segment.stop // chunkSize):
            counts[segment.chromosome][chunk] = 0
        # Normalizing counts
        if args.verbose:
            print('Normalising matches')
        normalizedCounts[direction] = segmentstatistics.normalizeCounts(counts)
        if args.verbose:
            print('Dumping normalized counts')
        countsFileName = '{}/ncounts.l{:02d}-{}.txt'.format(outputFolder, fragmentSize, direction)
        segmentutils.dumpNCounts(countsFileName, normalizedCounts[direction], chunkSize)

        # VI. GROUPING BY CYTOBANDS
        if args.cytomap is not None:
            if args.verbose:
                print('Grouping counts by cytomap regions')
            cytomapCounts[direction] = segmentstatistics.chunksToCytomap(cytomap,
                                                                         normalizedCounts[direction], chunkSize)
            cytoCountsFileName = '{}/cytocounts.l{:02d}-{}.txt'.format(outputFolder, fragmentSize, direction)
            segmentutils.dumpCytoCouns(cytoCountsFileName, cytomap, cytomapCounts[direction])

    # Outputting merged counts
    if args.mergedirections:
        nCountsMergedFileName = '{}/ncounts.l{:02d}-merged.txt'.format(outputFolder, fragmentSize)
        mergedCounts = {
            chromosome: (normalizedCounts['dir'][chromosome] + normalizedCounts['rev'][chromosome]) / 2
            for chromosome in normalizedCounts['dir'].keys()}
        segmentutils.dumpNCounts(nCountsMergedFileName, mergedCounts, chunkSize)
        if args.cytomap is not None:
            cytoCountsMergedFileName = '{}/cytocounts.l{:02d}-merged.txt'.format(outputFolder, fragmentSize)
            mergedCytoCounts = (cytomapCounts['dir'] + cytomapCounts['rev']) / 2
            segmentutils.dumpCytoCouns(cytoCountsMergedFileName, cytomap, mergedCytoCounts)
