#!/usr/bin/env python3

import argparse
import sys
import os
from os.path import join, abspath, curdir
from collections import Counter

# Adding local path to import program modules
sys.path.append(abspath(join(curdir, "segmentanalysis")))
from segmentanalysis import utils
from segmentanalysis import fragmentssearch, fragmentstatistics, genomeinterval
from segmentanalysis.typedef import (
    DnaSequence,
    Genome,
    FragmentsCounter,
)
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
    help="Provide file with set of patterns to exclude from calculations"
    + "Fragments containing at least one blacklisted pattern will be excluded from calculations",
)
parser.add_argument(
    "-i",
    "--include",
    type=str,
    help="Provide file with set of patterns to `protect` from blacklisting"
    + "Fragment containing at least one `include` pattern will not be blacklisted"
    + "If --blacklist argument is not set, this option will be ignored",
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
genome: Genome = utils.readGenome(utils.openMaybeGzipped(args.fastaFileName))

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

# II. ITERATION THROUGH ALL FRAGMENT SIZES AND DIRECTIONS
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

        # IV. Blacklisting fragments
        if args.blacklist is not None:
            fragmentsBlackSet: set[DnaSequence] = set(utils.readList(args.blacklist))
            fragmentsIncludeSet: set[DnaSequence] = (
                set(utils.readList(args.include)) if args.include is not None else set()
            )
            oldFragmentsCount = len(fragments)
            fragments = [
                fr
                for fr in fragments
                if (not utils.isContainSubFragments(fr, fragmentsBlackSet))
                or utils.isContainSubFragments(fr, fragmentsIncludeSet)
            ]
            if args.verbose:
                print(f"Balcklisted {oldFragmentsCount - len(fragments)} fragments")

        if len(fragments) == 0:
            print(
                f"No fragments selected for region {args.segment}, size {fragmentSize}, direction {direction}"
            )
            continue

        # V. SEARCH OF MATCHING FRAGMENTS
        fragmentsCounter: FragmentsCounter = Counter(fragments)

        chrFragmentsPositions = fragmentssearch.searchFragmentsGenome(
            genome, fragmentsCounter, args.verbose
        )
        if not args.nodump and fragmentSize >= args.mindumpsize:
            if args.verbose:
                print("Dumping fragments to text file")
            fragmentsFileName = f"{outputFolder}/fragments.l{fragmentSize:02d}-{segment.start}-{direction}.txt"
            fragmentssearch.dumpFragmentsToFile(
                fragmentsFileName, fragmentsCounter, chrFragmentsPositions
            )

        # VI. CONVERTING FRAGMENT POSITIONS TO COUNTS BY CHUNKS

        counts = fragmentstatistics.fragmentsToCounts(
            genome, chrFragmentsPositions, fragmentsCounter, chunkSize
        )

        countsFileName = f"{outputFolder}/rawcounts.l{fragmentSize:02d}-{segment.start}-{direction}.txt"

        utils.dumpCounts(countsFileName, counts, chunkSize)

        if segmentGenome is genome:  # Removing counts inside fragment
            if args.verbose:
                chunkStart, chunkStop = (
                    segment.start // chunkSize,
                    segment.stop // chunkSize,
                )
                print(
                    f"Setting to zero counts for chunks {segment.chromosome}:{chunkStart}-{chunkStop}"
                )
            fragmentstatistics.excludeIntervalFromCounts(counts, segment, chunkSize)
        rawCounts[direction] = counts

        # Normalizing counts
        if args.verbose:
            print("Normalising matches")
        normalizedCounts[direction] = fragmentstatistics.normalizeCounts(counts)
        if args.verbose:
            print("Dumping normalized counts")
        countsFileName = f"{outputFolder}/ncounts.l{fragmentSize:02d}-{segment.start}-{direction}.txt"
        utils.dumpCounts(
            countsFileName, normalizedCounts[direction], chunkSize, countsFormat="12.6f"
        )

        # VI. GROUPING BY CYTOBANDS
        if args.cytomap is not None:
            if args.verbose:
                print("Grouping counts by cytomap regions")
            cytomapCounts[direction] = fragmentstatistics.countsToCytomap(
                cytomap, normalizedCounts[direction], chunkSize
            )
            cytoCountsFileName = f"{outputFolder}/cytocounts.l{fragmentSize:02d}-{segment.start}-{direction}.txt"
            utils.dumpCytoCouns(cytoCountsFileName, cytomap, cytomapCounts[direction])

    # Outputting merged counts

    if args.verbose:
        print("Dumping merged counts")
    countsFileName = (
        f"{outputFolder}/rawcounts.l{fragmentSize:02d}-{segment.start}-merged.txt"
    )

    mergedRawCounts = fragmentstatistics.mergeCounts(rawCounts["dir"], rawCounts["rev"])

    utils.dumpCounts(countsFileName, mergedRawCounts, chunkSize, countsFormat="12.6f")

    nCountsMergedFileName = (
        f"{outputFolder}/ncounts.l{fragmentSize:02d}-{segment.start}-merged.txt"
    )

    mergednCounts = fragmentstatistics.mergeCounts(
        normalizedCounts["dir"], normalizedCounts["rev"]
    )
    utils.dumpCounts(
        nCountsMergedFileName, mergednCounts, chunkSize, countsFormat="12.6f"
    )

    if args.cytomap is not None:
        cytoCountsMergedFileName = (
            f"{outputFolder}/cytocounts.l{fragmentSize:02d}-{segment.start}-merged.txt"
        )
        mergedCytoCounts = (cytomapCounts["dir"] + cytomapCounts["rev"]) / 2
        utils.dumpCytoCouns(cytoCountsMergedFileName, cytomap, mergedCytoCounts)
