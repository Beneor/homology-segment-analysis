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


# 0. ANALYZING INPUT PARAMETERS

parser = argparse.ArgumentParser(description=programDescription, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("chromosomeseq", type=str,
                    help="TXT file containing chromosome sequence (one string, small letters, no name, no spaces)")
parser.add_argument("segment", type=str,
                    help="segment of chromosome to analyze (16113516:16900779) in BED-file notation (starting from 0, end is not included)")
parser.add_argument("-s", "--fragmentsizes", type=str, default='5,10,15,20,25,30',
                    help='Set of fragment sizes to search')
parser.add_argument("-d", "--fragmentdensity", type=int, default=25,
                    help='Average frequency of fragments in letters Default 25 meand that in average each 25-th letter will be start of fragment')

parser.add_argument("-c", "--chunk", type=int, default=10000,
                    help='Chunk size to divide chromosome')
parser.add_argument("-i", "--iterations", type=int, default=10,
                    help='Number of iterations for fragment splitting')
args = parser.parse_args()

# Processign command-line arguments
start, stop = [int(coord) for coord in args.segment.split(':')]
fragmentSizes = [int(size) for size in args.fragmentsizes.split(',')]

# I. INPUT FILES PREPARATION

# 1. Selection of chromosome segment 0 (to be fragmented)
with open(args.chromosomeseq) as chrSeqFile:
    chromosome=''.join( [line.strip() for line in chrSeqFile.readlines()] ).lower()
    segment = chromosome[start:stop]
#    f2 = open(r"Segment.txt", "w")  # File containing segment sequence
#    f2.write(segment)
#    f2.close()

# 2. Preparation of complementary segment
segmentRevComp = revcomp(segment)
#f = open(r"Segment.txt", "r")
#f2 = open(r"c-Segment.txt", "w")


# 3. Chromosome segmentation by chunks

# Fragments files mask
fileSuffix = '-x-fr.txt'
filesMask = '*' + fileSuffix
resultFiles = glob.glob(filesMask)

begin = 0  # the firts part start

# Sequence segmentation
fragments = open(r"x-small.txt", "r")  # Chromosome sequence to be segmented

for line in fragments:
    x = str(line)
    totallength = len(line)  # the length of x (in bases)
    seglength = args.chunk # the length of chunk (in bases)
    number = int(totallength / seglength)  # the number of parts
    for i in range(1, number + 1):
        outputFileName = str(i) + fileSuffix
        outputFile = open(outputFileName, "w")  # File with fragments
        end = begin + seglength  # chromosome part end
        seq = x[begin:end]  # chromosome part sequence
        outputFile.write(seq)  # random fragment sequence output
        begin += seglength  # the new part start
        outputFile.close()
fragments.close()

# II. SEARCH OF MATCHING FRAGMENTS

# For the direct chromosome DNA chain:

for i in range(iterations):  # 10 iterations

    for fragmentSize in fragmentSizes:
        pass


    match = open(r"5-dir.txt", "a+")  # 5 b matching fragments numbers
    match.write("Iteration ")
    match.write(str(i))
    match.write("\n")
    match.close()
    exec (open("./Search scripts/5-dir-restricted.py").read())

for i in range(1, 11):  # 10 iterations
    match = open(r"10-dir.txt", "a+")  # 15 b matching fragments numbers
    match.write("Iteration ")
    match.write(str(i))
    match.write("\n")
    match.close()
    exec (open("./Search scripts/10-dir-restricted.py").read())

for i in range(1, 11):  # 10 iterations
    match = open(r"15-dir.txt", "a+")  # 15 b matching fragments numbers
    match.write("Iteration ")
    match.write(str(i))
    match.write("\n")
    match.close()
    fragments = open(r"15-fragments-d.txt", "a+")  # 15 b matching fragments sequences and positions
    fragments.write("Iteration ")
    fragments.write(str(i))
    fragments.write("\n")
    fragments.close()
    exec (open("./Search scripts/15-dir.py").read())
    fragments = open(r"15-fragments-d.txt", "a+")
    fragments.write("\n")
    fragments.close()

for i in range(1, 11):  # 10 iterations
    match = open(r"20-dir.txt", "a+")  # 15 b matching fragments numbers
    match.write("Iteration ")
    match.write(str(i))
    match.write("\n")
    match.close()
    fragments = open(r"20-fragments-d.txt", "a+")  # 15 b matching fragments sequences and positions
    fragments.write("Iteration ")
    fragments.write(str(i))
    fragments.write("\n")
    fragments.close()
    exec (open("./Search scripts/20-dir.py").read())
    fragments = open(r"20-fragments-d.txt", "a+")
    fragments.write("\n")
    fragments.close()

for i in range(1, 11):  # 10 iterations
    match = open(r"25-dir.txt", "a+")  # 15 b matching fragments numbers
    match.write("Iteration ")
    match.write(str(i))
    match.write("\n")
    match.close()
    fragments = open(r"25-fragments-d.txt", "a+")  # 15 b matching fragments sequences and positions
    fragments.write("Iteration ")
    fragments.write(str(i))
    fragments.write("\n")
    fragments.close()
    exec (open("./Search scripts/25-dir.py").read())
    fragments = open(r"25-fragments-d.txt", "a+")
    fragments.write("\n")
    fragments.close()

for i in range(1, 11):  # 10 iterations
    match = open(r"30-dir.txt", "a+")  # 15 b matching fragments numbers
    match.write("Iteration ")
    match.write(str(i))
    match.write("\n")
    match.close()
    fragments = open(r"30-fragments-d.txt", "a+")  # 15 b matching fragments sequences and positions
    fragments.write("Iteration ")
    fragments.write(str(i))
    fragments.write("\n")
    fragments.close()
    exec (open("./Search scripts/30-dir.py").read())
    fragments = open(r"30-fragments-d.txt", "a+")
    fragments.write("\n")
    fragments.close()

# For the reverse chromosome DNA chain:

for i in range(1, 11):  # 10 iterations
    match = open(r"5-rev.txt", "a+")  # 15 b matching fragments numbers
    match.write("Iteration ")
    match.write(str(i))
    match.write("\n")
    match.close()
    exec (open("./Search scripts/5-rev-restricted.py").read())

for i in range(1, 11):  # 10 iterations
    match = open(r"10-rev.txt", "a+")  # 15 b matching fragments numbers
    match.write("Iteration ")
    match.write(str(i))
    match.write("\n")
    match.close()
    exec (open("./Search scripts/10-rev-restricted.py").read())

for i in range(1, 11):  # 10 iterations
    match = open(r"15-rev.txt", "a+")  # 15 b matching fragments numbers
    match.write("Iteration ")
    match.write(str(i))
    match.write("\n")
    match.close()
    fragments = open(r"15-fragments-r.txt", "a+")  # 15 b matching fragments sequences and positions
    fragments.write("Iteration ")
    fragments.write(str(i))
    fragments.write("\n")
    fragments.close()
    exec (open("./Search scripts/15-rev.py").read())
    fragments = open(r"15-fragments-r.txt", "a+")
    fragments.write("\n")
    fragments.close()

for i in range(1, 11):  # 10 iterations
    match = open(r"20-rev.txt", "a+")  # 15 b matching fragments numbers
    match.write("Iteration ")
    match.write(str(i))
    match.write("\n")
    match.close()
    fragments = open(r"20-fragments-r.txt", "a+")  # 15 b matching fragments sequences and positions
    fragments.write("Iteration ")
    fragments.write(str(i))
    fragments.write("\n")
    fragments.close()
    exec (open("./Search scripts/20-rev.py").read())
    fragments = open(r"20-fragments-r.txt", "a+")
    fragments.write("\n")
    fragments.close()

for i in range(1, 11):  # 10 iterations
    match = open(r"25-rev.txt", "a+")  # 15 b matching fragments numbers
    match.write("Iteration ")
    match.write(str(i))
    match.write("\n")
    match.close()
    fragments = open(r"25-fragments-r.txt", "a+")  # 15 b matching fragments sequences and positions
    fragments.write("Iteration ")
    fragments.write(str(i))
    fragments.write("\n")
    fragments.close()
    exec (open("./Search scripts/25-rev.py").read())
    fragments = open(r"25-fragments-r.txt", "a+")
    fragments.write("\n")
    fragments.close()

for i in range(1, 11):  # 10 iterations
    match = open(r"30-rev.txt", "a+")  # 15 b matching fragments numbers
    match.write("Iteration ")
    match.write(str(i))
    match.write("\n")
    match.close()
    fragments = open(r"30-fragments-r.txt", "a+")  # 15 b matching fragments sequences and positions
    fragments.write("Iteration ")
    fragments.write(str(i))
    fragments.write("\n")
    fragments.close()
    exec (open("./Search scripts/30-rev.py").read())
    fragments = open(r"30-fragments-r.txt", "a+")
    fragments.write("\n")
    fragments.close()

# III. OUTPUT FILES PREPARATION FOR ANALYSIS (SORTING)
exec (open(
    "./Sorting scripts/Sorting-initial-fragments-dir.py").read())  # Sorts all random (matching and dismatching) fragments according to their positions in segment
# (direct DNA chain). Input file:*-dir.all.txt, output file: *-dir-all.sort.txt
exec (open("./Sorting scripts/Sorting-initial-fragments-rev.py").read())  # Does the same for reverse DNA chain
exec (open(
    "./Sorting scripts/Sorting-matching-fragments-dir.py").read())  # Sorts matching fragments according to their positions in segment and chromosome (direct DNA chain)
# Input file: *-d.txt, output file: *-d-list.txt *-d-count.txt
exec (open("./Sorting scripts/Sorting-matching-fragments-rev.py").read())  # Does the same for reverse DNA chain
exec (open(
    "./Sorting scripts/Sorting-chromosome-segments-dir.py").read())  # Sorts numbers of matching fragments according to chromosome segment numbers
# and calculates average number normalized to the total fragments number (direct DNA chain)
exec (open("./Sorting scripts/Sorting-chromosome-segments-rev.py").read())  # Does the same for reverse DNA chain

# IV. ANALYSIS
# Table with the normalized average numbers of matching segments
exec (open(
    "./Analysis scripts/Normalization-table-dir.py").read())  # Generates table of normalized fragment frequencies of each length for each chromosome segment (direct DNA strain)
exec (open("./Analysis scripts/Normalization-table-rev.py").read())  # Does the same for reverse DNA chain
exec (open(
    "./Analysis scripts/Normalization-column-dir.py").read())  # Generates normalized fragment frequencies for each chromosomal disc and calculates Pearson correlation between fragment frequencies and ectopic contacts frequencies (direct DNA strain)
exec (open("./Analysis scripts/Normalization-column-rev.py").read())  # Does the same for reverse DNA chain
exec (open(
    "./Analysis scripts/Column-combinations.py").read())  # Combines columns with ectopic contacts and fragment frequencies for both DNA chains
exec (open(
    "./Analysis scripts/Normalization-column-both.py").read())  # Calculates Pearson correlation between fragment frequencies and ectopic contacts frequencies (both DNA strain)
