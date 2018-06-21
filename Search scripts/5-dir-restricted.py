#!/usr/bin/env python3

# Modules to be used
import random
import glob

# Chromosomal parts files mask
fileSuffix='-fr.txt'
filesMask='*'+fileSuffix
resultFiles = glob.glob(filesMask)

# I. Fragmentation of sequence 1
sequence = open(r"Segment.txt", "r") # Sequence 1 (to be fragmented)
fragmentedSeq = open(r"Words.txt", "w") # Sequence 1 fragments

for line in sequence:
        x=str(line)
        x1 = int(len(line)/25)# The length of lines in fragmented sequence / density of fragment choice: 1/25
        for i in range (1,x1):
                begin = random.randrange(len(line)) # random fragment beginning
                end= begin+5 #random fragment end
                seq=x[begin:end] # random fragment sequence
                print(seq)
                fragmentedSeq.write(seq) # random fragment sequence output
                fragmentedSeq.write("\n")
fragmentedSeq.close()
sequence.close()

# Seeking fragments in chromosomal parts

fragmentSeq = open(r"Words.txt", "r") # Sequence 1 fragments
fragmentNumber = open(r"5-dir.txt", "a+") # File containing the total count of matching fragments
print("%-12s%-12s%-17s"%("Matching N","Total N","File name"))
fragmentNumber.write("%-12s%-12s%-17s"%("Matching N","Total N","File name"))
fragmentNumber.write("\n")

for fileName in resultFiles:
        chromosomePart=open(fileName, "r")  #Chromosomal parts files
        for line in chromosomePart:
                x = str(line)
                matching = 0 # Number of matching fragments for chromosomal part
                total = 0 # Total number of fragments
                fragmentSeq.seek(0)
                for line in fragmentSeq:
                        a=str(line)
                        if len(a.strip())!=5: continue # checking the length of random fragment
                        matching+=x.count(a.strip()) #Calculating the total number of matching fragments
                        total+=1
        print("%-12s%-12s%-17s"%(matching,total,chromosomePart.name))
        fragmentNumber.write("%-12s%-12s%-17s"%(matching,total,chromosomePart.name)) #Writing the number of matching and total fragments for the chromosome segment
        fragmentNumber.write("\n")
        chromosomePart.close()
fragmentSeq.close()
fragmentNumber.close()

