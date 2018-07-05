# This module contains all code ralated to search of small fragments in chromosome chunks
import random
import re


class fragmentInfo:
    count = 0


def chooseFragments(segment, fragmentLength, fragmentDensity):
    """
    Generates list of random fragments form provided sequence
    :param segment: string sequence of segment
    :param fragmentLength: length of fragment
    :param fragmentDensity: density of fragments (see command-line parameters help
    :return: list of generated random fragments
    """
    fragments = []
    fragmentCount = len(segment) // fragmentDensity
    for i in range(fragmentCount):
        begin = random.randrange(len(segment) - fragmentLength)  # random fragment beginning
        end = begin + fragmentLength  # random fragment end
        fragment = segment[begin:end]  # random fragment sequence
        fragments.append(fragment)
    return fragments


def searchFragment(sequence, fragment):
    fragmentRe = re.compile(fragment)
    positions = [m.start() for m in fragmentRe.finditer(sequence)]
    return positions

def countFragments(sequence, fragments, verbose=False):
    nFragment = 0
    counts = {}
    nFragmentTotal = len(fragments)
    for fragment in fragments:
        fragmentCoordinates = searchFragment(sequence, fragment)
        counts[fragment] = len(fragmentCoordinates)
        nFragment += 1
        if verbose and nFragment % 100 == 0:
            print('Processed {}/{} fragments'.format(nFragment,nFragmentTotal))
    return counts