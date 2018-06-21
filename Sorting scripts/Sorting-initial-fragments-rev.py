# This script sorts matching fragments of segment(reverse sequence)
# according to the fragments position in segment

import glob

a = "New iteration"

#1. Sorting segment 1 fragments

fileSuffix='-rev-all.txt'
newfileSuffix='-rev-all-sort.txt'
filesMask='*'+fileSuffix
resultFiles = glob.glob(filesMask)

def sortfragments(line): #Sorting lines according to the files numbers in strings
        e = str(line)
        return int(e[int(e.find(";")+1):int(e.rfind(" "))].strip())#The first segment position

for fileName in resultFiles:
        matchingFragments=open(fileName, "r")  # Files with fragments - initial
        sortedMatchingFragments=open(fileName.rstrip(fileSuffix)+newfileSuffix,"w") # Sorted file
        print("%-12s"%(matchingFragments.name))
        for line in matchingFragments:
                x = str(line)
                bands = list()
                if a in x:
                        sortedMatchingFragments.write("New iteration")
                        for line in matchingFragments:
                                z=str(line)
                                if z.strip()=='': break
                                bands.append(line) #Creating list containing all strings
                        bands.sort(key=sortfragments) #Sorting values in the list
                        for i in bands:
                                sortedMatchingFragments.write(i) #Writing valuse of the list as separate strings
        sortedMatchingFragments.close()
        matchingFragments.close()   

