import io
from dataclasses import dataclass
from typing import TextIO
from segmentanalysis.typedef import Genome, DnaSequence


@dataclass(repr=False)
class GenomeInterval:
    """
    This class represents information about one genome interval in BED notation:
    chromosome, start, stop, optional interval name and numeric score
    """

    chromosome: DnaSequence
    start: int
    stop: int
    name: str = ""
    score: int = 0

    def __repr__(self) -> str:
        return (
            f"{self.chromosome}\t{self.start}\t{self.stop}\t{self.name}\t{self.score}"
        )


GenomeFeatures = list[GenomeInterval]


def strToBed(line: str, separator: str = "\t"):
    """
    Converts text line to BED interval
    :param line: string containing bed interval
    :param separator: separator used to split values in line
    :return: GenomeInterval instance containing readed interval
    """
    fields: list[str] = line.strip().split(separator)
    chromosome = fields[0]
    start, stop = [int(coord) for coord in fields[1:3]]
    intervalID = fields[3] if len(fields) > 3 else ""
    value = int(fields[4]) if len(fields) > 4 else 0

    return GenomeInterval(chromosome, start, stop, intervalID, value)


def readBedFile(bedFile: TextIO) -> GenomeFeatures:
    """
    Reads Information from file in BED notation. only 5 BED columns are supported
    :param bedFile: file IO - BED file to read
    :return: GenomeFeatures - list of GenomeInterval objects
    """
    genomeFeatures = []
    for line in bedFile:
        if line.strip()[0] == "#":  # Skipping comment
            continue
        interval = strToBed(line)
        genomeFeatures.append(interval)
    return genomeFeatures


def writeBedFile(bedFile: TextIO, genomeFeatures: GenomeFeatures) -> None:
    """
    Reads Information from file in BED notation. only 5 BED columns are supported
    :param bedFile: file IO - BED file to write
    :param genomeFeatures: GenomeFeatures object to write
    """
    for interval in genomeFeatures:
        bedFile.write(str(interval) + "\n")


def segmentStrToGGenomeInterval(segmentStr: str) -> GenomeInterval:
    """
    Converts bed-like location string to genome interval. Start and stop values can be omitted
    :param segmentStr: text string containing genome location
            chromosome:start:stop in full notation
            chromosome - in short notation
    :return: GenomeInterval object containing location
    """
    # Determining do we have start and stop positions in input string
    chromosomePos = segmentStr.find(":")
    interval = (
        strToBed(segmentStr, separator=":")
        if chromosomePos > 0
        else GenomeInterval(segmentStr, 0, 1)
    )
    return interval


def validateInterval(genome: Genome, interval: GenomeInterval) -> None:
    """
    Validates interval to Genome: checks for chromosome validness and moves
    :param genome:
    :param interval:
    """
    if not interval.chromosome in genome.keys():
        raise KeyError(
            f"Unknown chromosome id for loaded genome: {interval.chromosome}"
        )
    if interval.start >= interval.stop:
        raise ValueError(
            f"Interval with negative length: {interval.chromosome}:{interval.start}:{interval.stop}"
        )
    maxLen = len(genome[interval.chromosome])
    if interval.stop > maxLen:
        interval.stop = maxLen
