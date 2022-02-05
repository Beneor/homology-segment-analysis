from typing import Union
import numpy.typing as npt
from numpy import int64, float64

ChromosomeID = str
DnaSequence = str
Genome = dict[ChromosomeID, DnaSequence]

FragmentPosition = tuple[DnaSequence, int]
FragmentsCounter = dict[DnaSequence, int]

FragmentsLocations = dict[DnaSequence, list[int]]
GenomeFragmentsLocations = dict[ChromosomeID, FragmentsLocations]

Counts = Union[npt.NDArray[int64], npt.NDArray[float64]]
GenomeCounts = dict[ChromosomeID, Counts]
