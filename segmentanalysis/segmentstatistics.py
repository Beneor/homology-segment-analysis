from collections import defaultdict,OrderedDict

def locationsToChunks(fragmentPositions, chunkSize):
    """
    Converts exact positions of found matches to chunk numbers
    :param fragmentPositions: exact positions of fragments
    :param chunkSize: length of chunk
    :return: dictionaly:
       - key: fragment sequence
       - value: list of cunk numbers - one number for each match
    """
    return {fragment : [position // chunkSize for position in positions] 
                for fragment,positions in fragmentPositions.items()} 

def normalizeCounts(fragmentsPositionsChunks):
    """
    Normalizes counts 
    :param fragmentsPositionsChunks: positions of matchesfragments
    :param chunkSize: length of chunk
    :return: dictionaly:
       - key: fragment sequence
       - value: list of cunk numbers - one number for each match
    """
    pass
    
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
         
    
