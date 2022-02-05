# This module contains set of small service functions
import re
import gzip
from typing import TextIO

from segmentanalysis.typedef import DnaSequence, Genome, Counts, GenomeCounts
from segmentanalysis.genomeinterval import GenomeFeatures


# === DNA sequence utilities ===


def complement(sequence: DnaSequence) -> DnaSequence:
    """
    Make compliment of input sequence
    :param sequence: DNA sequence to make compliment
    :return: Compliment DNA sequence
    """
    base_complement = {"a": "t", "c": "g", "g": "c", "t": "a", "n": "n"}
    return "".join([base_complement[base] for base in sequence])


def revcomp(sequence: DnaSequence) -> DnaSequence:
    """
    Make reverse-compliment of input sequence
    :param sequence: DNA sequence to reverse-compliment
    :return: Reverse-compliment DNA sequence
    """
    return complement(sequence)[::-1]


# === IO FASTA read/write


def openMaybeGzipped(fileName: str) -> TextIO:
    """
    Determines whether input file is gzip-compressed and returns file handler via gzip / text file correspondingly
    :param fileName : Name of file to open
    :return: file handler to read lines
    """
    return (
        gzip.open(fileName, "rt") if fileName.endswith(".gz") else open(fileName, "r")
    )


def readList(filename: str) -> list[str]:
    """
    Reads strings from file and returns it as a list
    filename: file name to read strings
    :return: list of strings read from file filename
    """
    return [line.strip() for line in open(filename).readlines()]


def readGenome(fastaFile: TextIO) -> Genome:
    """
    Reads genome from fasta file
    :param fastaFile : Input file to read sequences
    :return: dictionary {'chromosomeID':'ChromosomeSequence'}
    """
    genome = {}  # Dictionary of chromosomes
    chrId = ""
    chrSeq: list[str] = []
    for line in fastaFile:
        if line[0] == ">":  # Next fasta record
            if chrId != "":  # Dumping current fast record
                genome[chrId] = "".join(chrSeq).lower()
            chrIdMatch = re.search(r">([-_A-z0-9]+)\s|$", line)
            if chrIdMatch is not None:
                chrId = chrIdMatch.group(1)  # Extracting new chromosome ID
            else:
                raise RuntimeError("No valid chromosome ID for line: " + line)
            chrSeq = []
        else:
            chrSeq.append(line.strip())
    genome[chrId] = "".join(
        chrSeq
    ).lower()  # Adding to dictionary last readed chromosome
    return genome


def dumpCounts(
    countsFileName: str,
    nCounts: GenomeCounts,
    chunkSize: int,
    countsFormat: str = "10d",
    skipZeros: bool = False,
) -> None:
    """
    writes calculated counts of fragments to text file
    :param countsFileName: name of file to write
    :param nCounts: dictionary containing normalized counts data from normalizeCounts function
    :param chunkSize:  size of one chunk
    """
    countsFile = open(countsFileName, "w")
    countsFormatter = "{0}\t{1}\t{2}\t{0}chunk{3}\t{4:" + countsFormat + "}\n"
    for chromosome in nCounts.keys():
        for chunk, nCount in enumerate(nCounts[chromosome]):
            if skipZeros and nCount == 0:
                continue
            start, stop = chunk * chunkSize, (chunk + 1) * chunkSize - 1

            countString = countsFormatter.format(chromosome, start, stop, chunk, nCount)
            countsFile.write(countString)
    countsFile.close()


def dumpCytoCouns(
    cytoCountsFileName: str, cytomap: GenomeFeatures, cytoCounts: Counts
) -> None:
    """
    Writes cytomap counts to BED file : cytobans genome coordinates, band ID and counts
    :param cytoCountsFileName: name of file to write
    :param cytomap: list of GenomeInterval objects representing cytomap
    :param cytoCounts: Numpy array of counts in the same order, as in cytomap list
    :return:
    """
    cytoCountsFile = open(cytoCountsFileName, "w")
    for i, interval in enumerate(cytomap):
        cytoCountsStr = f"{interval.chromosome}\t{interval.start}\t{interval.stop}\t{interval.name}\t{cytoCounts[i]}\n"
        cytoCountsFile.write(cytoCountsStr)
