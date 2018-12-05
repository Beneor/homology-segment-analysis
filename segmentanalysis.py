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
from os.path import join, abspath, curdir

sys.path.append(abspath(join(curdir, 'segmentanalysis')))
from segmentanalysis import segmentutils
from segmentanalysis import segmentsearch,segmentstatistics

# 0. ANALYZING INPUT PARAMETERS

parser = argparse.ArgumentParser(description=programDescription, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("fastaFileName", type=str,
                    help="FASTA file (may be gzipped) containing genome sequence")
parser.add_argument("segment", type=str,
                    help="segment of chromosome to analyze (X:16113516:16900779) in BED-file notation (starting from 0, end is not included)")

parser.add_argument("-v", "--verbose", action='store_true',
                    help='Print additional information to stdout')
parser.add_argument("-s", "--fragmentsizes", type=str, default='5,10,15,20,25,30',
                    help='Set of fragment sizes to search')
parser.add_argument("-d", "--fragmentdensity", type=int, default=10,
                    help='Average frequency of fragments in letters. Default 25 means that in average each 25-th letter will be start of fragment')

parser.add_argument("-c", "--chunk", type=int, default=10,
                    help='Chunk size to divide chromosome, in kilobases')

parser.add_argument("--dump", action='store_true',default=True,
                    help='Save found fragments and positions to files (used only for fragments >= mindumpsize)')
parser.add_argument("--mindumpsize", type=int, default=10,
                    help='Minumum size of fragment to store in file')

args = parser.parse_args()

# I. INPUT FILES PREPARATION

# 0. Reading genome sequence from FASTA[.gz] file
if args.fastaFileName.endswith('.gz'):
    fastaFile = gzip.open(args.fastaFileName, 'rt')
else:
    fastaFile = open(args.fastaFileName,'r')

if args.verbose:
    print('Loading genome from FASTA file')
genome = segmentutils.readFasta(fastaFile)
fastaFile.close()

# 1. Selection of chromosome segment
segment = segmentutils.segmentStrToGGenomeInterval(args.segment, genome)

segmentSeqFor = genome[segment.chromosome][segment.start:segment.stop]
segmentSeqs = {'dir':segmentSeqFor, 'rev':segmentutils.revcomp(segmentSeqFor)}

# 2. Size of chunks and selected fragments
fragmentSizes = [int(size) for size in args.fragmentsizes.split(',')]
chunkSize = args.chunk * 1000

# II. SEARCH OF MATCHING FRAGMENTS

# Generating random fragments from selected region
for direction,segmentSeq in segmentSeqs.items():
    for fragmentSize in fragmentSizes:
        # This is main cycle since all sizes handled individually
        if args.verbose:
            print('Generating fragments in {} direction for size: {}'.format(direction, fragmentSize))
        fragments = segmentsearch.chooseFragments(segmentSeq, fragmentSize, args.fragmentdensity)

        # Searching for fragments in chromosome
        chrFragmentsPositions = segmentsearch.searchFragments(genome, fragments, args.verbose)
        if args.dump and fragmentSize >= args.mindumpsize:
            if args.verbose:
                print('Dumping fragments to text file')
            fragmentsFileName = '{}.fragments.l{:02d}-{}.txt'.format(args.fastaFileName, fragmentSize, direction)
            segmentutils.dumpFragmentsToFile(fragmentsFileName, chrFragmentsPositions)

        chrPositionsChunks = segmentstatistics.locationsToChunks(chrFragmentsPositions, chunkSize)
        counts = segmentstatistics.chunksToCounts(chrPositionsChunks, genome, chunkSize)
        # Removing counts inside fragment
        if args.verbose:
            print('Setting to zero counts for chunks',segment.chromosome,
                  ':' ,segment.start // chunkSize,'-',segment.stop // chunkSize)
        for chunk in range(segment.start // chunkSize, segment.stop // chunkSize):
            counts[segment.chromosome][chunk] = 0
        # Normalizing counts
        if args.verbose:
            print('Normalising matches')
        normalizedCounts = segmentstatistics.normalizeCounts(counts)
        if args.dump:
            if args.verbose:
                print('Dumping normalized counts')
            countsFileName = '{}.ncounts.l{:02d}-{}.txt'.format(args.fastaFileName, fragmentSize, direction)
            segmentutils.dumpCountsToFile(countsFileName, normalizedCounts, chunkSize)

        # Group chunks by chromosome regions
        
        #Pearson correlation

# IV. ANALYSIS
    
