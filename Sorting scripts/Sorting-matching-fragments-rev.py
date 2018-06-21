# This script arranges the fragments of different lengths matching to chromosomal parts
# according to chromosome part number,
# then arranges fragments according to their position in segment (*-d-fin.txt),
# then arranges fragments according to their total number in chromosome (*-d-count.txt)
#(for reverse DNA sequence)

def sortnumber(line): #Sorting fragments matching to chromosome segments according to their position in chromosome
    x=str(line)
    return int(x[0:x.find("-")])

def sortnumber2(line): #Sorting fragments according to their position in segment
        e = str(line)
        return int(e[int(e.rfind(": ")+1):int(e.find(" ", e.rfind(": ")+2))].strip())#The first segment position

def sortnumber3(line): #Sorting fragments according to their total number in chromosome
        x = str(line)
        return int(x[x.find(":")+1:].strip())

import glob
import os

a = "Iteration"
b = "fr.txt"

# I. Modification of fragment files: space insertions between strings

fileSuffix='-r.txt'
newfileSuffix='-r-new.txt'
filesMask='*'+fileSuffix
resultFiles = glob.glob(filesMask)

for fileName in resultFiles:
        initial=open(fileName,"r+")
        modified=open(fileName.rstrip(fileSuffix)+newfileSuffix,"w")
        for line in initial:
                x=str(line)
                if a in x:
                        modified.write("\n")
                        modified.write("New iteration\n")
                        modified.write(x)
                        modified.write("\n")
                elif b in x:
                        modified.write(x)
                else:
                        modified.write(x)
initial.close()
modified.close()


# II. Arrangement of matching fragments according to segment position in chromosome

fileSuffix2="-r-new.txt"
newfileSuffix2="-r-list.txt"
filesMask2="*"+fileSuffix2
resultFiles2=glob.glob(filesMask2)
for fileName2 in resultFiles2:
    fragment=open(fileName2,"r+")
    fragmentArr=open(fileName2.rstrip(fileSuffix2)+newfileSuffix2,"w+")
    for line in fragment:
        x = str(line)
        bands=list()
        for line in fragment:
            x = str(line)
            if b in x:
                string=x
                for line in fragment:
                    x=str(line)
                    if x.strip()=="":break
                    string+=x
                bands.append(string)
        bands.sort(key=sortnumber)
        for i in bands:
            a = int(i[0:i.find("-")])
            if 1197 < a < 1278:continue
            fragmentArr.write(str(i))
            fragmentArr.write("\n")          
    fragment.close()
    fragmentArr.close()
    
#III. Arrangement of matching fragments according to their position in segment

fileSuffix3='-r-list.txt'
newfileSuffix3='-r-fin.txt'
filesMask3='*'+fileSuffix3
resultFiles3 = glob.glob(filesMask3)

for fileName3 in resultFiles3:
        matchingFragments=open(fileName3, "r")  # Files with fragments - initial
        sortedMatchingFragments=open(fileName3.rstrip(fileSuffix3)+newfileSuffix3,"w") # Sorted file
        print("%-12s"%(matchingFragments.name))
        for line in matchingFragments:
                x = str(line)
                if b in x:
                        sortedMatchingFragments.write(x)
                        bands=list()
                        for line in matchingFragments:
                                z=str(line)
                                if z.strip()=='': break
                                bands.append(line) #Creating list containing all strings
                        bands.sort(key=sortnumber2) #Sorting values in the list
                        for i in bands:
                                sortedMatchingFragments.write(i) #Writing valuse of the list as separate strings
        sortedMatchingFragments.close()
        matchingFragments.close()

#IV. Arrangement of matching fragments according to their total number in chromosome

fileSuffix4='-r-fin.txt'
newfileSuffix4='-r-temp.txt'
newfileSuffix5='-r-count.txt'
filesMask4='*'+fileSuffix4
resultFiles4 = glob.glob(filesMask4)

for fileName4 in resultFiles4:
        fragment=open(fileName4,"r+")
        comparedFragment=open(fileName4,"r+")
        fragmentCount=open(fileName4.rstrip(fileSuffix4)+newfileSuffix4,"w+")
        for line in fragment:
                x=str(line)
                if b in x: continue
                word=x[0:x.find(":")]
                i=0
                comparedFragment.seek(0)
                for line in comparedFragment:
                        y=str(line)
                        if b in y: continue
                        comparedWord=y[0:y.find(":")]
                        wordNum=int(y[y.find(":")+1:y.find("Positions")].strip())
                        if word in comparedWord:
                                i+=wordNum
                fragmentCount.write(word), fragmentCount.write(":"),fragmentCount.write(str(i)),fragmentCount.write("\n")          
        fragmentCount.close()
        fragment.close()
        comparedFragment.close()
        fragment2=open(fileName4.rstrip(fileSuffix4)+newfileSuffix5,"w+")
        fragment3=open(fileName4.rstrip(fileSuffix4)+newfileSuffix4,"r")
        bands=list()
        for line in fragment3:
            bands.append(line)
        unique = set(bands) #To leave out multiplicity
        bands=list(unique)
        bands.sort(key=sortnumber3)
        for i in bands:
            fragment2.write(i)
        fragment3.close()
        fragment2.close()
        fragment.close()

# V. Deletion of temporal files

fileSuffix6='-r-new.txt' #Deletes -new files
filesMask6='*'+fileSuffix6
resultFiles6 = glob.glob(filesMask6)
for f in resultFiles6:
        os.remove(f)

fileSuffix7='-r-list.txt' #Deletes -list files
filesMask7='*'+fileSuffix7
resultFiles7 = glob.glob(filesMask7)
for f in resultFiles7:
        os.remove(f)
        
fileSuffix8='-r-temp.txt' #Deletes -temp files
filesMask8='*'+fileSuffix8
resultFiles8 = glob.glob(filesMask8)
for f in resultFiles8:
        os.remove(f)
        
del a,b
