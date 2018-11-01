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
parser.add_argument("-d", "--fragmentdensity", type=int, default=100,
                    help='Average frequency of fragments in letters. Default 25 means that in average each 25-th letter will be start of fragment')

parser.add_argument("-c", "--chunk", type=int, default=10,
                    help='Chunk size to divide chromosome, in kilobases')
parser.add_argument("-i", "--iterations", type=int, default=10,
                    help='Number of iterations for fragment splitting')

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

# 1. Selection of chromosome segment 0 (to be fragmented
chromosome = args.segment.split(':')[0]
start, stop = [int(coord) for coord in args.segment.split(':')[1:]]
if not chromosome in genome.keys():
    print("Unknown chromosome id for loaded genome: {}".format(chromosome))
    print("Valid chromosomes are: "+','.join(genome.keys()))
    exit(-1)
if (start >= stop) or (stop > len(genome[chromosome])):
    print("Segment coordinates {}:{} are incorrect or greater then chromosome {} size: {}".format(
        start,stop,chromosome, len(genome[chromosome])))
    exit(-2)
segment = genome[chromosome][start:stop]
segmentRevComp = segmentutils.revcomp(segment)

# 2. Size of chunks and selected fragments
fragmentSizes = [int(size) for size in args.fragmentsizes.split(',')]
chunkSize = args.chunk * 1000

# II. SEARCH OF MATCHING FRAGMENTS

# Generating random fragments from selected region
fragments = []
for fragmentSize in fragmentSizes:
    # This is main cycle since all sizes handled individually
    if args.verbose:
        print('Generating fragments for size: {}'.format(fragmentSize))
    fragments += segmentsearch.chooseFragments(segment, fragmentSize, args.fragmentdensity)
    fragments += segmentsearch.chooseFragments(segmentRevComp, fragmentSize, args.fragmentdensity)
    
    # Searching for fragments in chromosome 
    fragmentsPositions = segmentsearch.searchFragments(genome, fragments, args.verbose)
    
    #Dumping fragments
    if args.dump and fragmentSize >= args.mindumpsize:
        if args.verbose:
            print('Dumping fragments to text file')
        dumpFile = open('{}.fragments.l{:02d}.txt'.format(args.fastaFileName,fragmentSize),'w')
        segmentutils.dumpFragmentsToFile(dumpFile, fragmentsPositions) 
        dumpFile.close()

exit(0)

# Converting found positions for nomalization
if args.verbose:
    print('Converting found positions for normalization')
fragmentsPositionsChunks = segmentstatistics.locationsToChunks(fragmentPositions, chunkSize)
chunksDensity = segmentstatistics.countDensity(fragmentsPositionsChunks)
normalizedDensity = segmentstatistics.normalizeDensity(chunksDensity)
print(normalizedDensity)

# IV. ANALYSIS
    



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
