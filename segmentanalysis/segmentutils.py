# This module contains set of small service functions

import re
import gzip

class GenomeInterval:
    """
    This class represents information about one genome interval in BED notation:
    chromosome, start, stop and optional text ID
    """

    def __init__(self, chromosome, start, stop, ID=None, value=None):
        self.chromosome = chromosome
        self.start = start
        self.stop = stop
        self.ID = ID
        self.value = value


def strToBed(line, separator='\t'):
    """
    Converts text line to BED interval
    :param line: string containing bed interval
    :param separator: separator used to split values in line
    :return: GenomeInterval instance containing readed interval
    """
    fields = line.strip().split(separator)
    chromosome = fields[0]
    start, stop = [int(coord) for coord in fields[1:3]]
    intervalID = fields[3] if len(fields) > 3 else None
    value = fields[4] if len(fields) > 4 else None
    return GenomeInterval(chromosome, start, stop, intervalID, value)


def segmentStrToGGenomeInterval(segmentStr, genome):
    """
    Converts bed-like location string to genome interval. Start and stop values can be omitted
    :param segmentStr: text string contatining denome location
            chromosome:start:stop in full notation
            chromosome - in short notation
    :param genome: - dictionaly containing loaded genome, for which location should be taken
    :return: GenomeInterval object containing location
    """
    # Determining do we have start and stop positions in input string
    chromosomePos = segmentStr.find(':')
    chromosome = segmentStr[:chromosomePos] if chromosomePos > 0 else segmentStr  # Extracting chromosome
    # Checking for consitency with genome
    if not chromosome in genome.keys():
        print("Unknown chromosome id for loaded genome: {}".format(chromosome))
        print("Valid chromosomes are: " + ','.join(genome.keys()))
        exit(-1)
    # Creating interval from full location strigng or for chromosome
    interval = strToBed(segmentStr, separator=':') if chromosomePos > 0 else GenomeInterval(
        segmentStr, 0, len(genome[segmentStr]))

    if (interval.start >= interval.stop) or (interval.stop > len(genome[interval.chromosome])):
        print("Segment coordinates {}:{} are incorrect or greater then chromosome {} size: {}".format(
            interval.start, interval.stop, interval.chromosome, len(genome[interval.chromosome])))
        exit(-2)
    return interval


def complement(sequence):
    """
    Make compliment of input sequence
    :param sequence: DNA sequence to make compliment
    :return: Compliment DNA sequence
    """
    basecomplement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return ''.join([basecomplement[base] for base in sequence])


def revcomp(sequence):
    """
    Make reverse-compliment of input sequence
    :param sequence: DNA sequence to revcomp
    :return: Reverse-compliment DNA sequence
    """
    return complement(sequence)[::-1]


def openMaybeGzipped(fileName):
    """
    Determines theter input file is gzip-compressed and returns file handler via gzip / text file correspondingly
    :param fileName : Name of file to open
    :return: file handler to read lines
    """
    return gzip.open(fileName, 'rt') if fileName.endswith('.gz') else open(fileName, 'r')


def readFasta(fastaFile):
    """
    Reads genome from fasta file
    :param fastaFile : Input file to read sequences
    :return: dictionary {'chromosomeID':'ChromosomeSequence'}
    """
    genome = {}  # Dictionary of chromosomes
    chrId = ''
    chrSeq = []
    for line in fastaFile:
        if line[0] == '>':  # Next fasta record
            if chrId != '':  # Dumping current fast record
                genome[chrId] = ''.join(chrSeq).lower()
            chrId = re.search(r'>([-_A-z0-9]+)\s|$', line).group(1)  # Extracting new chromosome ID
            chrSeq = []
        else:
            chrSeq.append(line.strip())
    genome[chrId] = ''.join(chrSeq).lower()  # Adding to dictionary last readed chromosome
    return genome


def dumpFragmentsToFile(fragmentsFileName, fragmentsPositions):
    """
    Writes set of fragments to file
    dumpFile: file stream to write data
    fragmentsPositions: dictionary with fragments and positions generated by segmentsearch.searchFragments() function
    """
    fragmentsFile = open(fragmentsFileName, 'w')
    for chromosome in fragmentsPositions.keys():
        sortedFragments = sorted(fragmentsPositions[chromosome].keys())
        for fragment in sortedFragments:
            fragmentString = '{}\t{:35}\t {}\n'.format(chromosome, fragment,
                                                       ','.join(
                                                           [str(i) for i in fragmentsPositions[chromosome][fragment]]))
            fragmentsFile.write(fragmentString)
    fragmentsFile.close()


def dumpNCounts(countsFileName, nCounts, chunkSize):
    """
    writes calculated counts of fragments to text file
    :param countsFileName: name of file to write
    :param nCounts: dictionary containing normalized counts data from normalizeCounts function
    :param chunkSize:  size of one chunk
    """
    countsFile = open(countsFileName, 'w')
    for chromosome in nCounts.keys():
        for chunk, nCount in enumerate(nCounts[chromosome]):
            start, stop = chunk * chunkSize, (chunk + 1) * chunkSize - 1
            countString = '{0}\t{1}\t{2}\t{0}chunk{3}\t{4:10.5f}\n'.format(chromosome, start, stop, chunk, nCount)
            countsFile.write(countString)
    countsFile.close()


def readBedFile(BedFileName):
    """
    Reads Information from file in BED notation. only 5 BED columns are supported
    :param BedFileName: name of
    :return: cytomap - list of GenomeInterval objects
    """
    bedFile = open(BedFileName)
    genomeFeatures = []
    IDs = set()
    for nLine, line in enumerate(bedFile):
        if line.strip()[0] == '#':  # Skipping comment
            continue
        interval = strToBed(line)
        genomeFeatures.append(interval)
        if interval.ID in IDs:
            print("WARNING: Non-unique interval ID {} at line {}".format(interval.ID, nLine))
        IDs.add(interval.ID)

    return genomeFeatures


def dumpCytoCouns(cytoCountsFileName, cytomap, cytoCounts):
    """
    Writes cytomap counts to BED file : cytobans genome coordinates, band ID and counts
    :param cytoCountsFileName: name of file to write
    :param cytomap: list of GenomeInterval objects representing cytomap
    :param cytoCounts: Numpy array of counts in the same order, as in cytomap list
    :return:
    """
    cytoCountsFile = open(cytoCountsFileName, 'w')
    for i, interval in enumerate(cytomap):
        cytoCountsStr = '{}\t{}\t{}\t{}\t{}\n'.format(
            interval.chromosome, interval.start, interval.stop, interval.ID, cytoCounts[i])
        cytoCountsFile.write(cytoCountsStr)
