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
parser.add_argument("-d", "--fragmentdensity", type=int, default=25,
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
segmentChromosome = args.segment.split(':')[0]
segmentStart, segmentStop = [int(coord) for coord in args.segment.split(':')[1:]]
if not segmentChromosome in genome.keys():
    print("Unknown chromosome id for loaded genome: {}".format(chromosome))
    print("Valid chromosomes are: "+','.join(genome.keys()))
    exit(-1)
if (segmentStart >= segmentStop) or (segmentStop > len(genome[segmentChromosome])):
    print("Segment coordinates {}:{} are incorrect or greater then chromosome {} size: {}".format(
        segmentStart,segmentStop,segmentChromosome, len(genome[segmentChromosome])))
    exit(-2)
segment = genome[segmentChromosome][segmentStart:segmentStop]
segmentRevComp = segmentutils.revcomp(segment)

# 2. Size of chunks and selected fragments
fragmentSizes = [int(size) for size in args.fragmentsizes.split(',')]
chunkSize = args.chunk * 1000

# II. SEARCH OF MATCHING FRAGMENTS

# Generating random fragments from selected region
for fragmentSize in fragmentSizes:
    # This is main cycle since all sizes handled individually
    if args.verbose:
        print('Generating fragments for size: {}'.format(fragmentSize))
    forFragments = segmentsearch.chooseFragments(segment, fragmentSize, args.fragmentdensity)
    revFragments = segmentsearch.chooseFragments(segmentRevComp, fragmentSize, args.fragmentdensity)
    
    # Searching for fragments in chromosome 
    chrFragmentsPositionsFor = segmentsearch.searchFragments(genome, forFragments, args.verbose)
    chrFragmentsPositionsRev = segmentsearch.searchFragments(genome, revFragments, args.verbose)
    
    #Dumping fragments
    if args.dump and fragmentSize >= args.mindumpsize:
        if args.verbose:
            print('Dumping fragments to text file')
        dumpFile = open('{}.fragments.l{:02d}-for.txt'.format(args.fastaFileName,fragmentSize),'w')
        segmentutils.dumpFragmentsToFile(dumpFile, chrFragmentsPositionsFor) 
        dumpFile.close()
        dumpFile = open('{}.fragments.l{:02d}-rev.txt'.format(args.fastaFileName,fragmentSize),'w')
        segmentutils.dumpFragmentsToFile(dumpFile, chrFragmentsPositionsRev) 
        dumpFile.close()
    if args.verbose:
        print('Normalising matches')
    chrPositionsChunksFor = segmentstatistics.locationsToChunks(chrFragmentsPositionsFor, chunkSize)
    chrPositionsChunksRev = segmentstatistics.locationsToChunks(chrFragmentsPositionsRev, chunkSize)
    
    # Removent counts inside fragment
    for chunk in range(segmentStart//chunkSize, segmentStop//chunkSize):
        print("resetting for chunk ", chunk)
        chrPositionsChunksFor[segmentChromosome][chunk] = 0
        chrPositionsChunksRev[segmentChromosome][chunk] = 0
    print(chrPositionsChunksFor[segmentChromosome])
    
    nFragmentsFor = segmentstatistics.normalizeCounts(chrPositionsChunksFor, genome, chunkSize)
    print (nFragmentsFor)
    #nragmentsRev  = 
    exit(0)

# Converting found positions for nomalization
chunksDensity = segmentstatistics.countDensity(fragmentsPositionsChunks)
normalizedDensity = segmentstatistics.normalizeDensity(chunksDensity)
print(normalizedDensity)

# IV. ANALYSIS
    



# III. OUTPUT FILES PREPARATION FOR ANALYSIS (SORTING)

