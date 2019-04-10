import itertools

try:
    import numpy as np  # Fast and memory efficient array calculations
except ImportError as error:
    # Output expected ImportErrors.
    print(
        '''Cannot import Numpy python module for fast array operations
        http://www.numpy.org/
    
        Try to install it by running "pip install numpy"
        or refer to module documentation''')
    exit(0)


def locationsToChunks(chrFragmentsPositions, chunkSize):
    """
    Converts exact positions of found matches in each chromosome to chunk numbers
    :param fragmentPositions: dictionary from segmentsearch.searchFragments()
    :param chunkSize: length of chunk
    :return: dictionary:
       - key: chromosome ID
       - value: list of cunk numbers - one number for each match
    """
    chrPositionsChunks = {}
    for chromosome, fragmentPositions in chrFragmentsPositions.items():
        # Wiping information about fragment sequences and joining all positions together
        # Also converting to faster numpy array
        positions = np.fromiter(itertools.chain(*fragmentPositions.values()), dtype=np.int64)
        chrPositionsChunks[
            chromosome] = positions // chunkSize  # Dividing each element by chunk size using numpy syntax
    return chrPositionsChunks


def chunksToCounts(chrPositionsChunks, genome, chunkSize):
    """
    Counts how many matches there are in each chunk
    :param chrPositionsChunks: positions of matches fragments for each chromosome
    :param genome: dictionary of chromosomes and genome sequences
    :param chunkSize: size of one chunk in bases
    :return: dictionary:
       - key: chromosome
       - value: list of counts for each chunk
    """
    counts = {}
    for chromosome, chunkPositions in chrPositionsChunks.items():
        counts[chromosome] = np.bincount(
            chunkPositions)  # Counting occurence of each chunk number and storing it into new list
        nChunks = len(genome[chromosome]) // chunkSize + 1
        if nChunks > len(counts[
                             chromosome]):  # Extending length to cover whole chromosome (need if we have no matches in latest chunks)
            counts[chromosome].resize(nChunks)
    return counts


def excludeIntervalFromCounts(counts, interval, chunkSize):
    for chunk in range(interval.start // chunkSize, interval.stop // chunkSize):
        counts[interval.chromosome][chunk] = 0
        
def normalizeCounts(counts):
    """
    Normalizes counts by mean count value
    :param counts: dictionary:
       - key: chromosome
       - value: list of counts for each chunk
    :return: dictionary:
       - key: chromosome
       - value: list of normalized counts for each chunk
    """
    nCounts = {}
    for chromosome, countsList in counts.items():
        sumCounts = np.sum(countsList)
        nChunks = len(countsList)
        averageCounts = sumCounts / nChunks
        nCounts[chromosome] = countsList / averageCounts if averageCounts > 0 else np.zeros(nChunks)
    return nCounts


def countsToCytomap(cytomap, counts, chunkSize):
    """
    Sums couns in chunks by cytomap
    :param cytomap: cytomap readed from bed file -- list genome interval and cytoband ID
    :param nCounts: dictionary containing normalized counts data from normalizeCounts function
    :param chunkSize: size of one chunk to calculate start and stop chunks for cytoband
    :return: numpy array containing summarized counts in same order as in cytomap
    """
    cytoCounts = np.zeros(len(cytomap))
    for i, interval in enumerate(cytomap):
        if not interval.chromosome in counts.keys(): # genome doesn't contain this interval
            print("WARNING:  Genome doen't contain chromosome listened in cytomap: ", interval.chromosome)
            continue
        for chunk in range(interval.start // chunkSize, interval.stop // chunkSize):
            cytoCounts[i] += counts[interval.chromosome][chunk]

    return cytoCounts
