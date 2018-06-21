# This script arranges fragments of specific length matching to chromosomal parts
# according to their total number within chromosome (10 iterations)

import glob
import os

a = "New iteration"
b="fr.txt"

# Chromosomal analysis result files

def sortnumber(line): #Sorting lines according to the files numbers in strings
        x = str(line)
        return int(x[x.find(":")+1:].strip())


fragment=open(r"30-fragment-d-fin.txt","r+")
comparedFragment=open(r"30-fragment-d-fin.txt","r+")
fragmentCount=open(r"Count-2.txt","w+")
for line in fragment:
    x=str(line)
    if a in x: continue
    if b in x: continue
    word=x[0:x.find(":")]
    i=0
    comparedFragment.seek(0)
    for line in comparedFragment:
        y=str(line)
        if a in y: continue
        if b in y: continue
        comparedWord=y[0:y.find(":")]
        wordNum=int(y[y.find(":")+1:y.find("Positions")].strip())   
        if word in comparedWord:
            i+=wordNum
    fragmentCount.write(word), fragmentCount.write(":"),fragmentCount.write(str(i)),fragmentCount.write("\n")          
fragmentCount.close()
fragment.close()
comparedFragment.close()
fragment2=open(r"Count-fin-2.txt","w+")
fragment3=open(r"Count-2.txt","r+")
bands=list()
for line in fragment3:
    bands.append(line)
unique = set(bands) #To leave out multiplicity
bands=list(unique)
bands.sort(key=sortnumber)
for i in bands:
    fragment2.write(i)
fragment3.close()
os.remove(r"Count-30-r.txt")
fragment2.close()
                
                
