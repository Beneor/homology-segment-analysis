from collections import defaultdict,OrderedDict

def locationsToChunks(fragmentPositions, chunkSize):
    fragmentsPositionsChunks={}
    for fragment,positions in fragmentPositions.items():
        fragmentsPositionsChunks[fragment] = [position // chunkSize for position in positions]
    return fragmentsPositionsChunks
    
    
def countDensity(fragmentsPositionsChunks):
    '''
    Reverses dict - counts for each chunk density of fragments
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
