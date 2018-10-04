# This module contains set of small service functions

def complement(sequence):
    """
    Make compliment of input sequence
    :param sequence: DNA sequence to make compliment
    :return:
    """
    basecomplement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return ''.join([basecomplement[base] for base in sequence])

def revcomp(sequence):
    """
    Make reverse-compliment of input sequence
    :param sequence: DNA sequence to revcomp
    :return:
    """
    return complement(sequence)[::-1]

def dumpFragmentsToFile(dumpFile, fragmentsPositions, minFragmentSize=15):
    """
    Writes set of fragments to file
    dumpFile: file stream to write data
    fragmentsPositions: dictionary with fragments and positions 
    minFragmentSize: minumum fragment position to write 
    """
    sortedFragments=sorted(fragmentsPositions.keys())
    for fragment in sortedFragments:
        if len(fragment) >= minFragmentSize:
            fragmentString = '{:35}: {}\n'.format(fragment, ','.join([str(i) for i in fragmentsPositions[fragment]]))
            dumpFile.write(fragmentString)
