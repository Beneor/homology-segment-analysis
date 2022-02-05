import itertools
from typing import Any
from segmentanalysis.typedef import (
    DnaSequence,
    Genome,
    FragmentsLocations,
    GenomeFragmentsLocations,
    FragmentsCounter,
    Counts,
    GenomeCounts,
)
from segmentanalysis.genomeinterval import GenomeInterval, GenomeFeatures

try:
    import numpy as np  # Fast and memory efficient array calculations
    import numpy.typing as npt
except ImportError as error:
    # Output expected ImportErrors.
    print(
        """Cannot import Numpy python module for fast array operations
        http://www.numpy.org/

        Try to install it by running "pip install numpy"
        or refer to module documentation"""
    )
    exit(0)


def fragmentsToCounts(
    genome: Genome,
    genomeFragmentsLocations: GenomeFragmentsLocations,
    fragmentsCounter: FragmentsCounter,
    chunkSize: int,
) -> GenomeCounts:
    """
    Converts exact positions of found matches in each chromosome to counts in one cunks
    :param genome:
    :param genomeFragmentsLocations: dictionary from segmentsearch.searchFragments()
    :param fragmentsCounter:
    :param chunkSize: length of chunk
    :return: dictionary:
       - key: chromosome ID
       - value: list of cunk numbers - one number for each match
    """
    genomeCounts: GenomeCounts = {}
    for chromosome, fragmentPositions in genomeFragmentsLocations.items():
        # Calculating expected chromosome size in chunks
        nChunks = len(genome[chromosome]) // chunkSize + 1
        genomeCounts[chromosome] = np.zeros(nChunks, dtype=np.int64)
        for fragment, positions in fragmentPositions.items():
            # Converting fragment exact coordinates to chunk coordinates (using fat numpy operations)
            chunkPositions = np.array(positions, dtype=np.int64) // chunkSize
            # Converting positions to counts
            counts = np.bincount(chunkPositions)
            # Updating counts to fragment multiplicity in original segment
            counts = counts * fragmentsCounter[fragment]
            # Extending length to cover whole chromosome (need if we have no matches in latest chunks)
            if nChunks > len(counts):
                counts.resize(nChunks)
            genomeCounts[chromosome] += counts
    return genomeCounts


def normalizeCounts(counts: GenomeCounts) -> GenomeCounts:
    """
    Normalizes counts by mean count value
    :param counts: dictionary:
       - key: chromosome
       - value: list of counts for each chunk
    :return: dictionary:
       - key: chromosome
       - value: list of normalized counts for each chunk
    """
    nCounts: GenomeCounts = {}
    for chromosome, countsList in counts.items():
        sumCounts = np.sum(countsList)
        nChunks = len(countsList)
        averageCounts = sumCounts / nChunks
        nCounts[chromosome] = (
            countsList / averageCounts
            if averageCounts > 0
            else np.zeros(nChunks, dtype=np.float64)
        )
    return nCounts


def mergeCounts(dirCounts: GenomeCounts, revCounts: GenomeCounts) -> GenomeCounts:
    """
    Calculates average for dir and rev counts on chromosome
    :param dirCounts: counts in forward direction
    :param revCounts: counts in reverse direction
    :return average counts for both directions
    """
    mergedCounts: GenomeCounts = {
        chromosome: (dirCounts[chromosome] + revCounts[chromosome]) / 2
        for chromosome in dirCounts.keys()
    }
    return mergedCounts


def countsToCytomap(
    cytomap: GenomeFeatures, counts: GenomeCounts, chunkSize: int
) -> Counts:
    """
    Sums counts in chunks by cytomap
    :param cytomap: cytomap readed from bed file -- list genome interval and cytoband ID
    :param counts: dictionary containing normalized counts data from normalizeCounts function
    :param chunkSize: size of one chunk to calculate start and stop chunks for cytoband
    :return: numpy array containing summarized counts in same order as in cytomap
    """
    cytoCounts = np.zeros(len(cytomap))
    for i, interval in enumerate(cytomap):
        if (
            not interval.chromosome in counts.keys()
        ):  # genome doesn't contain this interval
            print(
                "WARNING:  Genome doesn't contain chromosome listened in cytomap: ",
                interval.chromosome,
            )
            continue
        for chunk in range(interval.start // chunkSize, interval.stop // chunkSize):
            cytoCounts[i] += counts[interval.chromosome][chunk]

    return cytoCounts
