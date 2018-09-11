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
import glob
import sys
from os.path import join, abspath, curdir

sys.path.append(abspath(join(curdir, 'segmentanalysis')))
from segmentanalysis.segmentutils import revcomp
from segmentanalysis import segmentsearch,segmentstatistics

# 0. ANALYZING INPUT PARAMETERS

parser = argparse.ArgumentParser(description=programDescription, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("chromosomeseq", type=str,
                    help="TXT file containing chromosome sequence (one string, small letters, no name, no spaces)")
parser.add_argument("segment", type=str,
                    help="segment of chromosome to analyze (16113516:16900779) in BED-file notation (starting from 0, end is not included)")
parser.add_argument("-v", "--verbose", action='store_true',
                    help='Print additional information to stdout')
parser.add_argument("-s", "--fragmentsizes", type=str, default='5,10,15,20,25,30',
                    help='Set of fragment sizes to search')
parser.add_argument("-d", "--fragmentdensity", type=int, default=100,
                    help='Average frequency of fragments in letters. Default 25 means that in average each 25-th letter will be start of fragment')

parser.add_argument("-c", "--chunk", type=int, default=10,
                    help='Chunk size to divide chromosome, in kilobases')
parser.add_argument("-i", "--iterations", type=int, default=10,
                    help='Number of iterations for fragment splitting')

parser.add_argument("--dump", action='store_true',
                    help='Save found fragments and positions to files (used only for fragments >= mindumpsize)')
parser.add_argument("--mindumpsize", type=int, default=10,
                    help='Minumum size of fragment to store in file')

args = parser.parse_args()

# Processign command-line arguments
start, stop = [int(coord) for coord in args.segment.split(':')]
fragmentSizes = [int(size) for size in args.fragmentsizes.split(',')]
chunkSize = args.chunk * 1000

# I. INPUT FILES PREPARATION

# 1. Selection of chromosome segment 0 (to be fragmented)
with open(args.chromosomeseq) as chrSeqFile:
    chromosome = ''.join( [line.strip() for line in chrSeqFile.readlines() ]).lower()
    segment = chromosome[start:stop]

# 2. Preparation of complementary segment
segmentRevComp = revcomp(segment)

# II. SEARCH OF MATCHING FRAGMENTS
for iteration in range(args.iterations):  # iterations of segment splitting
    if args.verbose:
        print('Starting split iteration: ', iteration+1)
    # Generating random fragments from selected region
    fragments = []
    if args.verbose:
        print('Generating fragments for sizes: '+args.fragmentsizes)
    for fragmentSize in fragmentSizes: 
        fragments += segmentsearch.chooseFragments(segment, fragmentSize, args.fragmentdensity)
        fragments += segmentsearch.chooseFragments(segmentRevComp, fragmentSize, args.fragmentdensity)
    # Searching for fragments in chromosome 
    fragmentPositions = segmentsearch.searchFragments(chromosome, fragments, args.verbose)
    
    # Converting found positions for nomalization
    if args.verbose:
        print('Converting found positions for normalization')
    fragmentsPositionsChunks = segmentstatistics.locationsToChunks(fragmentPositions, chunkSize)
    chunksDensity = segmentstatistics.countDensity(fragmentsPositionsChunks)
    print(chunksDensity)

# IV. ANALYSIS
    

exit(0)

# III. OUTPUT FILES PREPARATION FOR ANALYSIS (SORTING)
'''
exec(open(
    "./Sorting scripts/Sorting-initial-fragments-dir.py").read())  # Sorts all random (matching and dismatching) fragments according to their positions in segment
# (direct DNA chain). Input file:*-dir.all.txt, output file: *-dir-all.sort.txt

exec(open("./Sorting scripts/Sorting-initial-fragments-rev.py").read())  # Does the same for reverse DNA chain
exec(open(
    "./Sorting scripts/Sorting-matching-fragments-dir.py").read())  # Sorts matching fragments according to their positions in segment and chromosome (direct DNA chain)
# Input file: *-d.txt, output file: *-d-list.txt *-d-count.txt
exec(open("./Sorting scripts/Sorting-matching-fragments-rev.py").read())  # Does the same for reverse DNA chain
exec(open(
    "./Sorting scripts/Sorting-chromosome-segments-dir.py").read())  # Sorts numbers of matching fragments according to chromosome segment numbers
# and calculates average number normalized to the total fragments number (direct DNA chain)
exec(open("./Sorting scripts/Sorting-chromosome-segments-rev.py").read())  # Does the same for reverse DNA chain
'''

# Table with the normalized average numbers of matching segments
exec(open(
    "./Analysis scripts/Normalization-table-dir.py").read())  # Generates table of normalized fragment frequencies of each length for each chromosome segment (direct DNA strain)
exec(open("./Analysis scripts/Normalization-table-rev.py").read())  # Does the same for reverse DNA chain
exec(open(
    "./Analysis scripts/Normalization-column-dir.py").read())  # Generates normalized fragment frequencies for each chromosomal disc and calculates Pearson correlation between fragment frequencies and ectopic contacts frequencies (direct DNA strain)
exec(open("./Analysis scripts/Normalization-column-rev.py").read())  # Does the same for reverse DNA chain
exec(open(
    "./Analysis scripts/Column-combinations.py").read())  # Combines columns with ectopic contacts and fragment frequencies for both DNA chains
exec(open(
    "./Analysis scripts/Normalization-column-both.py").read())  # Calculates Pearson correlation between fragment frequencies and ectopic contacts frequencies (both DNA strain)
