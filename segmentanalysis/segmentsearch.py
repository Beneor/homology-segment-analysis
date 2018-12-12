# This module contains all code related to searching of small fragments in chromosome chunks
import random
import re
from collections import defaultdict

try:
    import ahocorasick as aho_corasick # Fast and memory efficient library for exact or approximate multi-pattern string search
except ImportError as error:
    # Output expected ImportErrors.
    print(
    '''Cannot import Aho-Corasick python module for fast fragment search
    https://pypi.org/project/pyahocorasick/
    
    Try to install it by running "pip install pyahocorasick"
    or refer to module documentation
    For Windows Visual studio buildtools are required: https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2017''')
    exit(0)

outputInfoChunkLength = 1000000 # Chromosome length chunk to print fragment finding statistics

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

def searchFragments(genome, fragments, verbose=False):
    """
    Finds multiple fragments in chromosome sequence
    :param genome: dictionaly of strings -- {chromosomeID:chromosomeSequence}
    :param fragments: list of fragments to search
    :return: dictionary of found fragments positions:
        - key: chromosome ID
        - value: dictionary:
            - key: fragment sequence 
            - value: list of positions in chromosome
    """
    if verbose:
        print('Making Aho-Corasick automation for {} fragments'.format(len(fragments)))
    A = aho_corasick.Automaton()
    for idx, key in enumerate(fragments):
        A.add_word(key, (idx, key))
    A.make_automaton()
    
    if verbose:
        print('Starting Aho-Corasick search for {} fragments'.format(len(fragments)))
    chrFragmentsPositions = {}
    
    for chrId,sequence in genome.items(): 
        chrFragmentsPositions[chrId] = defaultdict(list)
        if verbose:
            print('Processing chromosome '+chrId)
            currentInfoChunk = 0
        for end_index, (insert_order, fragment) in A.iter(sequence):
            position = end_index - len(fragment) + 1
            chrFragmentsPositions[chrId][fragment].append(position) 
            if verbose and (position//outputInfoChunkLength) > currentInfoChunk:
                currentInfoChunk = position//outputInfoChunkLength
                print('Processed chromosome {}, {:8}/{:8} bases'.format(chrId, position,len(sequence)))
    if verbose:
        print('Aho-Corasick search finished')
    return chrFragmentsPositions
    
