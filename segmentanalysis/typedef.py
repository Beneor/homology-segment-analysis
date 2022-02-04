ChromosomeID = str
DnaSequence = str
Genome = dict[ChromosomeID, DnaSequence]
FragmentPosition = tuple[DnaSequence, int]
FragmentsCounter = dict[DnaSequence, int]

FragmentsLocations = dict[DnaSequence, list[int]]

GenomeFragmentsLocations = dict[ChromosomeID, FragmentsLocations]



