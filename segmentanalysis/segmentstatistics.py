import itertools
from collections import Counter
import numpy as np

#from collections import defaultdict,OrderedDict

def locationsToChunks(chrFragmentsPositions, chunkSize):
    """
    Converts exact positions of found matches in each chromosome to chunk numbers
    :param fragmentPositions: dictionary from segmentsearch.searchFragments()
    :param chunkSize: length of chunk
    :return: dictionaly:
       - key: chromosome ID
       - value: list of cunk numbers - one number for each match
    """
    chrPositionsChunks={}
    for chromosome, fragmentPositions in chrFragmentsPositions.items():
        # Wiping information about fragment sequences and joining all positions together
        positions = np.fromiter(itertools.chain(*fragmentPositions.values()),dtype=np.int64)
        chrPositionsChunks[chromosome] = positions//chunkSize
    return chrPositionsChunks 

def normalizeCounts(chrPositionsChunks, genome, chunkSize):
    """
    Normalizes counts 
    :param chrPositionsChunksFor: positions of matchesfragments for each chromosome
    :return: dictionaly:
       - key: chromosome
       - value: list of count for each chunk
    """
    nFragments = {}
    for chromosome, chunkPositions in chrPositionsChunks.items():
        ncounts = np.bincount(chunkPositions) # Counting ocuurence of each cunk nomber and storing it into new list
        nChunks = len(genome[chromosome])/chunkSize
        print(chromosome)
        print(ncounts)
        print("sum of counts", np.sum(ncounts))
        print("chunks: ",nChunks)
        averageCouns = np.sum(ncounts) / nChunks
        print("average counts for each chunk: ", averageCouns)
        nFragments[chromosome] = ncounts/averageCouns
        
    return nFragments

    
def countDensity(fragmentsPositionsChunks):
    '''
    "Inverses" counts dict - counts for each chunk density of fragments
    Result of this fucntion is following:
     dict:
         key: is a chunk number
         value: dict:
             key: fragment length
             value: fragments count
    '''
    chunksDensity= {}
    for fragment,chunks in fragmentsPositionsChunks.items():
        for chunk in chunks:
            if not chunk in chunksDensity.keys():
                chunksDensity[chunk]=defaultdict(int)
            chunksDensity[chunk][len(fragment)]+=1
    
    return OrderedDict(sorted(chunksDensity.items()))

def normalizeDensity(chunksDensity):
    '''
    '''
    for chunk,fragmentsCount in chunksDensity.items():
        pass
         
    
