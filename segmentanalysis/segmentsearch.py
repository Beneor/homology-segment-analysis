# This module contains all code ralated to search of small fragments in chromosome chunks
import random
import re
import ahocorasick as aho_corasick # Fast and memory efficient library for exact or approximate multi-pattern string search
from collections import defaultdict

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

def countFragments(sequence, fragments, verbose=False):
    nFragment = 0
    counts = defaultdict(list)
    if verbose:
        print('Making Aho-Corasick automation for {} fragments'.format(len(fragments)))
    A = aho_corasick.Automaton()
    for idx, key in enumerate(fragments):
        A.add_word(key, (idx, key))
    A.make_automaton()
    if verbose:
        print('Starting Aho-Corasick search for {} fragments'.format(len(fragments)))
    foundFragments = 0
    for end_index, (insert_order, fragment) in A.iter(sequence):
        position = end_index - len(fragment) + 1
        counts[fragment].append(position) 
        foundFragments +=1
        if verbose and foundFragments % 100000 == 0:
            print('Found {} fragments, search position: {}'.format(foundFragments,position))
    if verbose:
        print('Aho-Corasick search finished')
    return counts

