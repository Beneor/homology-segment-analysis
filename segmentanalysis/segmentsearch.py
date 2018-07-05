# This module contains all code ralated to search of small fragments in chromosome chunks
import random


def chooseFragments(segment, fragmentLength, fragmentDensity):
    """
    Generates list of random fragments form provided sequence
    :param segment: string sequence of segment
    :param fragmentLength: length of fragment
    :param fragmentDensity: density of fragments (see command-line parameters help
    :return: list of generated random fragments
    """
    fragments = []
    fragmentCount = len(segment) / fragmentDensity
    for i in range(fragmentCount):
        begin = random.randrange(len(segment) - fragmentLength)  # random fragment beginning
        end = begin + fragmentLength  # random fragment end
        fragment = segment[begin:end]  # random fragment sequence
        fragments.append(fragment)
        print(fragment)
    return fragments
