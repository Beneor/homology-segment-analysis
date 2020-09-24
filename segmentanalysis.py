#!/usr/bin/env python3

import argparse
import sys
import os
from os.path import join, abspath, curdir

# Adding local path to import program modules
sys.path.append(abspath(join(curdir, 'segmentanalysis')))
from segmentanalysis import segmentutils
from segmentanalysis import segmentsearch, segmentstatistics

programDescription = '''
This script produces a set of short (20-40 by default) fragments for the specified sequence segment,
splits genome chromosomes into a set of long (10 kb by default) chunks and
seeks for matching fragments within chunks.
Then it computes frequencies of matching fragments fir each chunk and normalizes it.
Also it calculates fragments matching frequency for cytobands, it cytoband information is provided.

See details in program description in README.md  
'''

# 0. ANALYZING INPUT PARAMETERS
parser = argparse.ArgumentParser(description=programDescription, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("fastaFileName", type=str,
                    help="FASTA file (may be gzipped) containing genome sequence")
parser.add_argument("segment", type=str,
                    help="segment of chromosome to analyze: file and location " +
                         "(:X:11982050:12772070). See README.md for details")
parser.add_argument("cytomap", type=str, nargs='?',
                    help="BED file containing cytobands")

parser.add_argument("-v", "--verbose", action='store_true',
                    help='Print additional information to stdout')
parser.add_argument("-s", "--fragmentsizes", type=str, default='20,25,30,35,40',
                    help='Set of fragment sizes to search')
parser.add_argument("-d", "--fragmentdensity", type=float,
                    help='Use random chosen fragments instead of all possible. ' +
                         'Sets average frequency of fragments in letters. ' +
                         'For example 5 means that each 5-th letter will be start of fragment')
parser.add_argument("-c", "--chunk", type=float, default=10.0,
                    help='Chunk size to divide chromosome, in kilobases')
parser.add_argument("-b", "--blacklist", type=str,
                    help='Provide file with set of fragments to exclude from calculations')
parser.add_argument("-i", "--include", type=str,
                    help='Provide file with set of fragments to include into calculations')
parser.add_argument("-f", "--force", action='store_true',
                    help='Write to existing output data folder. Can overwrite old data')
parser.add_argument("--nodump", action='store_true',
                    help='Do not save found fragments and positions to files')
parser.add_argument("--mindumpsize", type=int, default=20,
                    help='Minimum size of fragment to store exact locations in file')

args = parser.parse_args()

# I. INPUT FILES PREPARATION

# 0. Reading genome sequence from FASTA[.gz] file
#if args.fastaFileName.endswith('.gz'):
##    fastaFile = gzip.open(args.fastaFileName, 'rt')
#else:
#    fastaFile = open(args.fastaFileName, 'r')

if args.verbose:
    print('Loading genome from FASTA file')
genome = segmentutils.readFasta(segmentutils.openMaybeGzipped(args.fastaFileName))

# 1. Selection of chromosome segment
locationPos = args.segment.index(':')
# Determining do we have seconn file name in location, and reading seconf genome if nesessary
segmentGenome = segmentutils.readFasta(open(args.segment[:locationPos], 'r')) if locationPos != 0 else genome
# Extracting segment location
segment = segmentutils.segmentStrToGGenomeInterval(args.segment[locationPos + 1:], segmentGenome)
segmentSeqFor = segmentGenome[segment.chromosome][segment.start:segment.stop]
segmentSeqs = {'dir': segmentSeqFor, 'rev': segmentutils.revcomp(segmentSeqFor)}

# 2. Size of chunks and selected fragments
fragmentSizes = [int(size) for size in args.fragmentsizes.split(',')]
chunkSize = int(args.chunk * 1000)

# 3. Folder to save results
outputFolder = '{}.{}-{}-{}'.format(args.fastaFileName, segment.chromosome, segment.start, segment.stop)
if os.path.exists(outputFolder):
    if not args.force:
        print('Output folder {} already exists\nPlease delete and resubmit'.format(outputFolder))
        exit(0)
else:
    os.makedirs(outputFolder)

if args.cytomap is not None:
    cytomap = segmentutils.readBedFile(args.cytomap)

if args.blacklist is not None:
    fragmentsBlackList = set(segmentutils.readList(args.blacklist))
    
if args.include is not None:
    fragmentsIncludeList = set(segmentutils.readList(args.include))

# II. ITERATION THROUGH ALL FRAGMENT SIZES AND DIRECTIONS

# Generating random fragments from selected region
for fragmentSize in fragmentSizes:
    normalizedCounts = {}
    cytomapCounts = {}
    for direction, segmentSeq in segmentSeqs.items():
        # III. GENERATING FRAGMENTS
        if args.verbose:
            print('Generating fragments in {} direction for size: {}'.format(direction, fragmentSize))
        if args.fragmentdensity is not None:
            fragments = segmentsearch.chooseRandomFragments(segmentSeq, fragmentSize, args.fragmentdensity)
        else:
            fragments = segmentsearch.makeFragments(segmentSeq, fragmentSize)
        
        if args.include is not None:
            if args.blacklist is not None:
                oldFragmentsCount = len(fragments)
                fragments = [fr for fr in fragments if not (segmentutils.isContainSubFragments(fr, fragmentsBlackList) and not segmentutils.isContainSubFragments(fr, fragmentsIncludeList))]
                if args.verbose:
                    print('Balcklisted and included {} fragments'.format(oldFragmentsCount - len(fragments)))
        else:
            if args.blacklist is not None:
                oldFragmentsCount = len(fragments)
                fragments = [fr for fr in fragments if not segmentutils.isContainSubFragments(fr, fragmentsBlackList)]
                if args.verbose:
                    print('Balcklisted and included {} fragments'.format(oldFragmentsCount - len(fragments)))

        # IV. SEARCH OF MATCHING FRAGMENTS
        chrFragmentsPositions = segmentsearch.searchFragments(genome, fragments, args.verbose)
        if not args.nodump and fragmentSize >= args.mindumpsize:
            if args.verbose:
                print('Dumping fragments to text file')
            fragmentsFileName = '{}/fragments.l{:02d}-{}-{}.txt'.format(outputFolder, fragmentSize, segment.start, direction)
            segmentutils.dumpFragmentsToFile(fragmentsFileName, chrFragmentsPositions)
        # v. CONVEERTING TO CHUNKS AND NORMALIZING
        chrPositionsChunks = segmentstatistics.locationsToChunks(chrFragmentsPositions, chunkSize)
        counts = segmentstatistics.chunksToCounts(chrPositionsChunks, genome, chunkSize)
        if segmentGenome is genome:  # Removing counts inside fragment
            if args.verbose:
                print('Setting to zero counts for chunks {}:{}-{}'.format(segment.chromosome,
                      segment.start // chunkSize, segment.stop // chunkSize))
            segmentstatistics.excludeIntervalFromCounts(counts, segment, chunkSize)
        # Normalizing counts
        if args.verbose:
            print('Normalising matches')
        normalizedCounts[direction] = segmentstatistics.normalizeCounts(counts)
        if args.verbose:
            print('Dumping normalized counts')
        countsFileName = '{}/ncounts.l{:02d}-{}-{}.txt'.format(outputFolder, fragmentSize, segment.start, direction,)
        segmentutils.dumpNCounts(countsFileName, normalizedCounts[direction], chunkSize)

        # VI. GROUPING BY CYTOBANDS
        if args.cytomap is not None:
            if args.verbose:
                print('Grouping counts by cytomap regions')
            cytomapCounts[direction] = segmentstatistics.countsToCytomap(cytomap,
                                                                         normalizedCounts[direction], chunkSize)
            cytoCountsFileName = '{}/cytocounts.l{:02d}-{}-{}.txt'.format(outputFolder, fragmentSize, segment.start, direction)
            segmentutils.dumpCytoCouns(cytoCountsFileName, cytomap, cytomapCounts[direction])

    # Outputting merged counts
    nCountsMergedFileName = '{}/ncounts.l{:02d}-{}-merged.txt'.format(outputFolder,fragmentSize,segment.start)
    #Corrected filename!!!!
    mergedCounts = {
        chromosome: (normalizedCounts['dir'][chromosome] + normalizedCounts['rev'][chromosome]) / 2
        for chromosome in normalizedCounts['dir'].keys()}
    segmentutils.dumpNCounts(nCountsMergedFileName, mergedCounts, chunkSize)
    if args.cytomap is not None:
        cytoCountsMergedFileName = '{}/cytocounts.l{:02d}-{}-merged.txt'.format(outputFolder, fragmentSize,segment.start)
        mergedCytoCounts = (cytomapCounts['dir'] + cytomapCounts['rev']) / 2
        segmentutils.dumpCytoCouns(cytoCountsMergedFileName, cytomap, mergedCytoCounts)
