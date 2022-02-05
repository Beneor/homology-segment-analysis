# This module contains all code related to searching of small fragments in chromosome chunks
import random
from typing import Iterator
from collections import defaultdict
from segmentanalysis.typedef import (
    DnaSequence,
    Genome,
    FragmentPosition,
    GenomeFragmentsLocations,
    FragmentsCounter,
)
from segmentanalysis.consts import outputInfoChunkLength

# Fast and memory efficient library for exact or approximate multi-pattern string search
try:
    import ahocorasick as aho_corasick  # type: ignore
except ImportError as error:
    # Output expected ImportErrors.
    print(
        """Cannot import Aho-Corasick python module for fast fragment search
        https://pypi.org/project/pyahocorasick/

        Try to install it by running "pip3 install pyahocorasick"
        or refer to module documentation
        For Windows Visual studio buildtools are required: https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2017"""
    )
    exit(0)


def makeFragments(segment: DnaSequence, fragmentLength: int) -> list[DnaSequence]:
    """
    Generates set of all possible unique fragments from provided sequence
    :param segment: segment string to split
    :param fragmentLength: size of fragment to generate
    :return: list of generated fragments
    """
    # Collecting all fragments with potential duplicates
    fragments = [
        segment[start : start + fragmentLength]
        for start in range(len(segment) - fragmentLength)
    ]
    return fragments


def chooseRandomFragments(
    segment: DnaSequence, fragmentLength: int, fragmentDensity: float
) -> list[DnaSequence]:
    """
    Generates set of unique random fragments from provided sequence
    :param segment: string sequence of segment
    :param fragmentLength: length of fragment
    :param fragmentDensity: density of fragments (see command-line parameters help
    :return: set of generated random fragments
    """
    fragmentCount = int(len(segment) / fragmentDensity)
    # Collecting all fragments with potential duplicates
    fragments = [""] * fragmentCount
    for i in range(fragmentCount):
        # random fragment beginning
        begin = random.randrange(len(segment) - fragmentLength)
        end = begin + fragmentLength  # random fragment end
        fragments[i] = segment[begin:end]  # random fragment sequence
    return fragments


def searchFragmentsChromosome(
    automation: aho_corasick.Automaton, chromosome: DnaSequence
) -> Iterator[FragmentPosition]:
    """
    Makes Aho-Corasick search for the provided chromosome sequence
    :param automation: Aho-Corasick automation
    :param chromosome: string of one chromosome
    :return: Iterator with positions of fragments
    """
    for end_index, (insert_order, fragment) in automation.iter(chromosome):
        position = end_index - len(fragment) + 1
        yield fragment, position


def searchFragmentsGenome(
    genome: Genome, fragmentsCounter: FragmentsCounter, verbose: bool = False
) -> GenomeFragmentsLocations:
    """
    Finds multiple fragments in chromosome sequence
    :param genome: dictionaly of strings -- {chromosomeID:chromosomeSequence}
    :param fragmentsCounter: Counter dict containing fragments with counts
    :param verbose: print log messages
    :return: dictionary of found fragments positions:
        - key: chromosome ID
        - value: dictionary:
            - key: fragment sequence
            - value: list of positions in chromosome
    """
    if verbose:
        print(f"Making Aho-Corasick automation for {len(fragmentsCounter)} fragments")
    A = aho_corasick.Automaton()
    for idx, key in enumerate(
        fragmentsCounter.keys()
    ):  # Adding to automation all unique fragments
        A.add_word(key, (idx, key))
    A.make_automaton()

    if verbose:
        print(f"Starting Aho-Corasick search for {len(fragmentsCounter)} fragments")
    chrFragmentsPositions: GenomeFragmentsLocations = {}

    for chrId, sequence in genome.items():
        chrFragmentsPositions[chrId] = defaultdict(list)
        if verbose:
            print(f"Processing chromosome {chrId}")
            currentInfoChunk = 0
        for fragment, position in searchFragmentsChromosome(A, sequence):
            chrFragmentsPositions[chrId][fragment].append(position)
            if verbose and (position // outputInfoChunkLength) > currentInfoChunk:
                currentInfoChunk = position // outputInfoChunkLength
                print(
                    f"Processed chromosome {chrId}, {position:8}/{len(sequence):8} bases"
                )
    if verbose:
        print("Aho-Corasick search finished")
    return chrFragmentsPositions


def dumpFragmentsToFile(
    fragmentsFileName: str,
    fragmentsCounter: FragmentsCounter,
    fragmentsPositions: GenomeFragmentsLocations,
) -> None:
    """
    Writes set of fragments to file
    dumpFile: file stream to write data
    fragmentsPositions: dictionary with fragments and positions generated by segmentsearch.searchFragments() function
    """
    fragmentsFile = open(fragmentsFileName, "w")
    for chromosome in fragmentsPositions.keys():
        sortedFragments = sorted(fragmentsPositions[chromosome].keys())
        for fragment in sortedFragments:
            positions = ",".join(
                [str(i) for i in fragmentsPositions[chromosome][fragment]]
            )
            fragmentString = (
                f"{chromosome}\t{fragment}\t{fragmentsCounter[fragment]}\t{positions}\n"
            )
            fragmentsFile.write(fragmentString)
    fragmentsFile.close()
