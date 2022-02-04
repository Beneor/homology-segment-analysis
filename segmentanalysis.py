#!/usr/bin/env python3

import argparse
import sys
import os
from os.path import join, abspath, curdir
from collections import Counter

# Adding local path to import program modules
sys.path.append(abspath(join(curdir, "segmentanalysis")))
from segmentanalysis import utils
from segmentanalysis import fragmentssearch, segmentstatistics, genomeinterval
from segmentanalysis.typedef import DnaSequence, Genome, FragmentsCounter
from segmentanalysis.genomeinterval import GenomeInterval

programDescription = """
This script produces a set of short (20-40 by default) fragments for the specified sequence segment,
splits genome chromosomes into a set of long (10 kb by default) chunks and
seeks for matching fragments within chunks.
Then it computes frequencies of matching fragments fir each chunk and normalizes it.
Also it calculates fragments matching frequency for cytobands, it cytoband information is provided.

See details in program description in README.md
"""

# 0. ANALYZING INPUT PARAMETERS
parser = argparse.ArgumentParser(
    description=programDescription, formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument(
    "fastaFileName",
    type=str,
    help="FASTA file (may be gzipped) containing genome sequence",
)
parser.add_argument(
    "segment",
    type=str,
    help="segment of chromosome to analyze: file and location "
    + "(:X:11982050:12772070). See README.md for details",
)
parser.add_argument(
    "cytomap", type=str, nargs="?", help="BED file containing cytobands"
)

parser.add_argument(
    "-v",
    "--verbose",
    action="store_true",
    help="Print additional information to stdout",
)
parser.add_argument(
    "-s",
    "--fragmentsizes",
    type=str,
    default="20,25,30,35,40",
    help="Set of fragment sizes to search",
)
parser.add_argument(
    "-d",
    "--fragmentdensity",
    type=float,
    help="Use random chosen fragments instead of all possible. "
    + "Sets average frequency of fragments in letters. "
    + "For example 5 means that each 5-th letter will be start of fragment",
)
parser.add_argument(
    "-c",
    "--chunk",
    type=float,
    default=10.0,
    help="Chunk size to divide chromosome, in kilobases",
)
parser.add_argument(
    "-b",
    "--blacklist",
    type=str,
    help="Provide file with set of fragments to exclude from calculations",
)
parser.add_argument(
    "-i",
    "--include",
    type=str,
    help="Provide file with set of fragments to include into calculations",
)
parser.add_argument(
    "-f",
    "--force",
    action="store_true",
    help="Write to existing output data folder. Can overwrite old data",
)
parser.add_argument(
    "--nodump",
    action="store_true",
    help="Do not save found fragments and positions to files",
)
parser.add_argument(
    "--mindumpsize",
    type=int,
    default=20,
    help="Minimum size of fragment to store exact locations in file",
)

args = parser.parse_args()

# I. INPUT FILES PREPARATION

if args.verbose:
    print("Loading genome from FASTA file")
genome: Genome = utils.readGenome(
    utils.openMaybeGzipped(args.fastaFileName)
)

# 1. Selection of chromosome segment
locationPos = args.segment.index(":")
# Determining do we have second file name in location, and reading second genome if necessary
segmentGenome: Genome = (
    utils.readGenome(open(args.segment[:locationPos], "r"))
    if locationPos != 0
    else genome
)

segment: GenomeInterval = genomeinterval.segmentStrToGGenomeInterval(
    args.segment[locationPos + 1 :]
)
try:
    genomeinterval.validateInterval(segmentGenome, segment)
except Exception as err:
    print(err)
    exit(-1)

# Extracting segment location
segmentSeqFor: str = segmentGenome[segment.chromosome][segment.start : segment.stop]
segmentSeqs: dict[str, str] = {
    "dir": segmentSeqFor,
    "rev": utils.revcomp(segmentSeqFor),
}

# 2. Size of chunks and selected fragments
fragmentSizes: list[int] = [int(size) for size in args.fragmentsizes.split(",")]
chunkSize = int(args.chunk * 1000)

# 3. Folder to save results
outputFolder = "{}.{}-{}-{}".format(
    args.fastaFileName, segment.chromosome, segment.start, segment.stop
)
if os.path.exists(outputFolder):
    if not args.force:
        print(
            f"Output folder {outputFolder} already exists\nPlease delete and resubmit"
        )
        exit(0)
else:
    os.makedirs(outputFolder)

if args.cytomap is not None:
    with utils.openMaybeGzipped(args.cytomap) as cytomapFile:
        cytomap = genomeinterval.readBedFile(cytomapFile)

if args.blacklist is not None:
    fragmentsBlackList: set[DnaSequence] = set(utils.readList(args.blacklist))

if args.include is not None:
    fragmentsIncludeList: set[DnaSequence] = set(utils.readList(args.include))

# II. ITERATION THROUGH ALL FRAGMENT SIZES AND DIRECTIONS

# Generating random fragments from selected region
for fragmentSize in fragmentSizes:
    rawCounts = {}
    normalizedCounts = {}
    cytomapCounts = {}
    for direction, segmentSeq in segmentSeqs.items():
        # III. GENERATING FRAGMENTS
        if args.verbose:
            print(
                f"Generating fragments in {direction} direction for size: {fragmentSize}"
            )
        if args.fragmentdensity is not None:
            fragments = fragmentssearch.chooseRandomFragments(
                segmentSeq, fragmentSize, args.fragmentdensity
            )
        else:
            fragments = fragmentssearch.makeFragments(segmentSeq, fragmentSize)

        if args.include is not None:
            if args.blacklist is not None:
                oldFragmentsCount = len(fragments)
                fragments = [
                    fr
                    for fr in fragments
                    if not (
                        utils.isContainSubFragments(fr, fragmentsBlackList)
                        and not utils.isContainSubFragments(
                            fr, fragmentsIncludeList
                        )
                    )
                ]
                if args.verbose:
                    print(
                        f"Balcklisted and included {oldFragmentsCount - len(fragments)} fragments"
                    )
        else:
            if args.blacklist is not None:
                oldFragmentsCount = len(fragments)
                fragments = [
                    fr
                    for fr in fragments
                    if not utils.isContainSubFragments(fr, fragmentsBlackList)
                ]
                if args.verbose:
                    print(
                        f"Balcklisted and included {oldFragmentsCount - len(fragments)} fragments"
                    )

        if len(fragments) == 0:
            print(
                f"No fragments selected for region {args.segment}, size {fragmentSize}, direction {direction}"
            )
            continue

        # IV. SEARCH OF MATCHING FRAGMENTS

        fragmentsCounter: FragmentsCounter = Counter(fragments)

        chrFragmentsPositions = fragmentssearch.searchFragmentsGenome(
            genome, fragmentsCounter, args.verbose
        )
        if not args.nodump and fragmentSize >= args.mindumpsize:
            if args.verbose:
                print("Dumping fragments to text file")
            fragmentsFileName = f"{outputFolder}/fragments.l{fragmentSize:02d}-{segment.start}-{direction}.txt"
            fragmentssearch.dumpFragmentsToFile(fragmentsFileName, fragmentsCounter, chrFragmentsPositions)

        # v. CONVEERTING TO CHUNKS AND NORMALIZING
        chrPositionsChunks = segmentstatistics.locationsToChunks(
            chrFragmentsPositions, chunkSize
        )
        counts = segmentstatistics.chunksToCounts(chrPositionsChunks, genome, chunkSize)

        if segmentGenome is genome:  # Removing counts inside fragment
            if args.verbose:
                chunkStart, chunkStop = (
                    segment.start // chunkSize,
                    segment.stop // chunkSize,
                )
                print(
                    f"Setting to zero counts for chunks {segment.chromosome}:{chunkStart}-{chunkStop}"
                )
            segmentstatistics.excludeIntervalFromCounts(counts, segment, chunkSize)
        rawCounts[direction] = counts

        # Normalizing counts
        if args.verbose:
            print("Normalising matches")
        normalizedCounts[direction] = segmentstatistics.normalizeCounts(counts)
        if args.verbose:
            print("Dumping normalized counts")
        countsFileName = f"{outputFolder}/ncounts.l{fragmentSize:02d}-{segment.start}-{direction}.txt"
        utils.dumpCounts(
            countsFileName, normalizedCounts[direction], chunkSize, countsFormat="10.4f"
        )

        # VI. GROUPING BY CYTOBANDS
        if args.cytomap is not None:
            if args.verbose:
                print("Grouping counts by cytomap regions")
            cytomapCounts[direction] = segmentstatistics.countsToCytomap(
                cytomap, normalizedCounts[direction], chunkSize
            )
            cytoCountsFileName = f"{outputFolder}/cytocounts.l{fragmentSize:02d}-{segment.start}-{direction}.txt"
            utils.dumpCytoCouns(
                cytoCountsFileName, cytomap, cytomapCounts[direction]
            )

    # Outputting merged counts

    if args.verbose:
        print("Dumping merged raw counts")
    countsFileName = (
        f"{outputFolder}/rawcounts.l{fragmentSize:02d}-{segment.start}-merged.txt"
    )

    mergedRawCounts = {
        chromosome: (rawCounts["dir"][chromosome] + rawCounts["rev"][chromosome])
        for chromosome in rawCounts["dir"].keys()
    }

    utils.dumpCounts(
        countsFileName, mergedRawCounts, chunkSize, countsFormat="10d", skipZeros=True
    )

    nCountsMergedFileName = (
        f"{outputFolder}/ncounts.l{fragmentSize:02d}-{segment.start}-merged.txt"
    )
    # Corrected filename!!!!
    mergedCounts = {
        chromosome: (
            normalizedCounts["dir"][chromosome] + normalizedCounts["rev"][chromosome]
        )
        / 2
        for chromosome in normalizedCounts["dir"].keys()
    }
    utils.dumpCounts(
        nCountsMergedFileName, mergedCounts, chunkSize, countsFormat="10.4f"
    )

    if args.cytomap is not None:
        cytoCountsMergedFileName = (
            f"{outputFolder}/cytocounts.l{fragmentSize:02d}-{segment.start}-merged.txt"
        )
        mergedCytoCounts = (cytomapCounts["dir"] + cytomapCounts["rev"]) / 2
        utils.dumpCytoCouns(cytoCountsMergedFileName, cytomap, mergedCytoCounts)
