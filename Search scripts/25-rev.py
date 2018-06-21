#!/usr/bin/env python3

# Modules to be used
import random
import glob

# Chromosomal parts files mask
fileSuffix='-fr.txt'
filesMask='*'+fileSuffix
resultFiles = glob.glob(filesMask)

# I. Fragmentation of sequence 1
sequence = open(r"c-Segment.txt", "r") # Sequence 1 (to be fragmented)
fragmentedSeq = open(r"Words.txt", "w") # Sequence 1 fragments
allFragmentedSeq = open(r"25-words-rev-all.txt", "a+") #Sequence 1 fragments - all iterations
allFragmentedSeq.write(" \n")
allFragmentedSeq.write("New iteration")
allFragmentedSeq.write("\n")

for line in sequence:
        x=str(line)
        x1 = int(len(line)/25)# The length of lines in fragmented sequence / density of fragment choice: 1/25
        for i in range (1,x1):
                begin = random.randrange(len(line)) # random fragment beginning
                end= begin+25 #random fragment end
                seq=x[begin:end] # random fragment sequence
                print(seq)
                fragmentedSeq.write(seq) # random fragment sequence output
                fragmentedSeq.write("\n")
                allFragmentedSeq.write(seq) # random fragment sequence output
                allFragmentedSeq.write(";") # random fragment sequence output
                allFragmentedSeq.write(str(begin))
                allFragmentedSeq.write("\n")
fragmentedSeq.close()
sequence.close()

# Seeking fragments in chromosomal parts

fragmentSeq = open(r"Words.txt", "r") # Sequence 1 fragments
fragmentPosition = open(r"25-fragments-r.txt", "a+") # File containing fragment sequences, number and positions of matches
fragmentNumber = open(r"25-rev.txt", "a+") # File containing the total count of matching fragments
print("%-12s%-12s%-17s"%("Matching N","Total N","File name"))
fragmentNumber.write("%-12s%-12s%-17s"%("Matching N","Total N","File name"))
fragmentNumber.write("\n")

for fileName in resultFiles:
        chromosomePart=open(fileName, "r")  #Chromosomal parts files
        fragmentPosition.write("%-12s"%(chromosomePart.name))
        fragmentPosition.write("\n")
        for line in chromosomePart:
                x = str(line)
                matching = 0 # Number of matching fragments for chromosomal part
                total = 0 # Total number of fragments
                fragmentSeq.seek(0)
                for line in fragmentSeq:
                        a=str(line)
                        if len(a.strip())!=25: continue # checking the length of random fragment
                        total+=1 #Calculating the total number of fragments
                        if x.count(a.strip())!=0: #checking that fragment is within the analyzed sequence
                                fragmentPosition.write(a.strip()) #Writing fragment sequence
                                fragmentPosition.write(": ")
                                fragmentPosition.write(str(x.count(a.strip()))) #Writing the begin of the first matching fragment within the chromosome segment
                                fragmentPosition.write(" Positions: ")
                                start=0 #Begin of the scanned area of chromosome segment
                                i=1
                                while True: # output of fragment positions for multiple matches
                                        if i>x.count(a.strip()): 
                                                break                   
                                        start = x.find(a.strip(),start) #Finding the begin of the scanned area
                                        fragmentPosition.write(str(x.find(a.strip(),start))) #Writing the begin of the next matching fragment within the chromosome segment
                                        fragmentPosition.write(str(" "))
                                        start+=25 #Shifting the begin of the scanned area
                                        i+=1
                                fragmentPosition.write(str("\n"))
                        matching+=x.count(a.strip()) #Calculating the total number of matching fragments
                fragmentPosition.write(" \n")
        print("%-12s%-12s%-17s"%(matching,total,chromosomePart.name))
        fragmentNumber.write("%-12s%-12s%-17s"%(matching,total,chromosomePart.name)) #Writing the number of matching and total fragments for the chromosome segment
        fragmentNumber.write("\n")
        chromosomePart.close()
fragmentSeq.close()
fragmentPosition.close()
fragmentNumber.close()
