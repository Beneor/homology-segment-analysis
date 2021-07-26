# This module contains all code related to searching of small fragments in chromosome chunks
import random
from collections import defaultdict,Counter

warnFragmentCount = 500

try:
    import \
        ahocorasick as aho_corasick  # Fast and memory efficient library for exact or approximate multi-pattern string search
except ImportError as error:
    # Output expected ImportErrors.
    print(
        '''Cannot import Aho-Corasick python module for fast fragment search
        https://pypi.org/project/pyahocorasick/
        
        Try to install it by running "pip install pyahocorasick"
        or refer to module documentation
        For Windows Visual studio buildtools are required: https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2017''')
    exit(0)

outputInfoChunkLength = 1000000  # Chromosome length chunk to print fragment finding statistics


def makeFragments(segment, fragmentLength):
    """
    Generates set of all possible unique fragments from provided sequence
    :param segmentSeq: segment string to split
    :param fragmentSize: sife of fragment to generate
    :return: set of generated fragments
    """
    fragments = [] # Collecting all fragments with potential duplicates
    for start in range(len(segment) - fragmentLength):
        fragments.append(segment[start:start + fragmentLength])
    return fragments


def chooseRandomFragments(segment, fragmentLength, fragmentDensity):
    """
    Generates set of unique random fragments from provided sequence
    :param segment: string sequence of segment
    :param fragmentLength: length of fragment
    :param fragmentDensity: density of fragments (see command-line parameters help
    :return: set of generated random fragments
    """
    fragments = [] # Collecting all fragments with potential duplicates
    fragmentCount = int(len(segment) / fragmentDensity)
    for i in range(fragmentCount):
        begin = random.randrange(len(segment) - fragmentLength)  # random fragment beginning
        end = begin + fragmentLength  # random fragment end
        fragments.append(segment[begin:end])  # random fragment sequence
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
    fragmentsCounter = Counter(fragments) # Converting to counter since Aho-Corasic works only with unique fragments
    for idx, key in enumerate(fragmentsCounter.keys()): # Adding to automation all uqique fragments
        A.add_word(key, (idx, key))
    A.make_automaton()
    
    if verbose:
        print('Starting Aho-Corasick search for {} fragments'.format(len(fragments)))
    chrFragmentsPositions = {}
    for chrId, sequence in genome.items():
        chrFragmentsPositions[chrId] = defaultdict(list)
        if verbose:
            print('Processing chromosome ' + chrId)
            currentInfoChunk = 0
        for end_index, (insert_order, fragment) in A.iter(sequence):
            position = end_index - len(fragment) + 1
            chrFragmentsPositions[chrId][fragment].append(position)
            if verbose and (position // outputInfoChunkLength) > currentInfoChunk:
                currentInfoChunk = position // outputInfoChunkLength
                print('Processed chromosome {}, {:8}/{:8} bases'.format(chrId, position, len(sequence)))
        
        for fragment in chrFragmentsPositions[chrId]: # Checking for duplicated fragments
            if fragmentsCounter[fragment] > 1: # Duplicating positions list for each repeated fragment
                if fragmentsCounter[fragment] > warnFragmentCount:
                    print("! WARNING: Overpresented k-mer. May cause high memory consumption")
                    print("k-mer {}:{} counts".format(fragment, fragmentsCounter[fragment]))
                chrFragmentsPositions[chrId][fragment] = chrFragmentsPositions[chrId][fragment] * fragmentsCounter[fragment]
    if verbose:
        print('Aho-Corasick search finished')
    return chrFragmentsPositions
