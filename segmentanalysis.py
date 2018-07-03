#!/usr/bin/env python3
# This script produces a set of short (5 - 30 b) random segments for the specified chromosome part,
# particulates chromosome into a set of long (10 kb) parts,
# seeks random segments within chromosome parts,
# computes frequencies of occurence for segments of different length within chromosome parts,
# makes table for the normalized frequencies,
# and calculates Pearson correlation between fragment frequencies for chromosome discs and ectopic contact frequencies


# I. INPUT FILES PREPARATION

# 1. Selection of chromosome segment 0 (to be fragmented)
f = open(r"x-small.txt", "r", encoding="cp1251") #Chromosome sequence (one string, small letters, no name, no spaces)

for line in f:
        x = str(line)
        y=x[16113516:16900779] #14B-15B segment borders. !!! CHOOSE THE NUMBERS OF SEGMENT BEGIN (M-1) AND END (N)
        f2=open(r"Segment.txt", "w", encoding="cp1251") #File containing segment sequence
        f2.write(y)
        f2.close()
f.close()

# 2. Preparation of complementary segment
f = open(r"Segment.txt", "r")
f2 = open(r"c-Segment.txt","w")

for line in f:
        x = str(line)
        for i in range(0,len(x)):
                j = int(len(x)-i-1)
                s = x[j]
                if s == "a":
                        f2.write("t")
                elif s == "g":
                        f2.write("c")
                elif s == "c":
                        f2.write("g")
                elif s == "t":
                        f2.write("a")
                else:
                        f2.write("n")

# 3. Chromosome segmentation

# Modules to be used
import random
import glob

# Fragments files mask
fileSuffix='-x-fr.txt'
filesMask='*'+fileSuffix
resultFiles = glob.glob(filesMask)

length = 10 # the length of chromosome segments (in kilobases)
begin = 0 # the firts part start

# Sequence segmentation
fragments = open(r"x-small.txt", "r") #Chromosome sequence to be segmented

for line in fragments:
        x=str(line)
        totallength=len(line) #the length of x (in bases)
        seglength = int(length)*1000 #the length of part (in bases)
        number=int(totallength/seglength) #the number of parts
        for i in range(1,number+1):
                outputFileName=str(i)+fileSuffix
                outputFile=open(outputFileName, "w") # File with fragments
                end=begin+seglength #chromosome part end
                seq=x[begin:end] #chromosome part sequence
                outputFile.write(seq) #random fragment sequence output
                begin+=seglength #the new part start
                outputFile.close()
fragments.close()

# II. SEARCH OF MATCHING FRAGMENTS

#For the direct chromosome DNA chain:

for i in range (1,11): # 10 iterations
        match = open(r"5-dir.txt", "a+") #15 b matching fragments numbers
        match.write ("Iteration ")
        match.write (str(i))
        match.write ("\n")
        match.close()
        exec(open("./Search scripts/5-dir-restricted.py").read())

for i in range (1,11): # 10 iterations
        match = open(r"10-dir.txt", "a+") #15 b matching fragments numbers
        match.write ("Iteration ")
        match.write (str(i))
        match.write ("\n")
        match.close()
        exec(open("./Search scripts/10-dir-restricted.py").read())

for i in range (1,11): # 10 iterations
        match = open(r"15-dir.txt", "a+") #15 b matching fragments numbers
        match.write ("Iteration ")
        match.write (str(i))
        match.write ("\n")
        match.close()
        fragments = open(r"15-fragments-d.txt", "a+") #15 b matching fragments sequences and positions
        fragments.write ("Iteration ")
        fragments.write (str(i))
        fragments.write ("\n")
        fragments.close()
        exec(open("./Search scripts/15-dir.py").read())
        fragments = open(r"15-fragments-d.txt", "a+")
        fragments.write ("\n")
        fragments.close()

for i in range (1,11): # 10 iterations
        match = open(r"20-dir.txt", "a+") #15 b matching fragments numbers
        match.write ("Iteration ")
        match.write (str(i))
        match.write ("\n")
        match.close()
        fragments = open(r"20-fragments-d.txt", "a+") #15 b matching fragments sequences and positions
        fragments.write ("Iteration ")
        fragments.write (str(i))
        fragments.write ("\n")
        fragments.close()
        exec(open("./Search scripts/20-dir.py").read())
        fragments = open(r"20-fragments-d.txt", "a+")
        fragments.write ("\n")
        fragments.close()

for i in range (1,11): # 10 iterations
        match = open(r"25-dir.txt", "a+") #15 b matching fragments numbers
        match.write ("Iteration ")
        match.write (str(i))
        match.write ("\n")
        match.close()
        fragments = open(r"25-fragments-d.txt", "a+") #15 b matching fragments sequences and positions
        fragments.write ("Iteration ")
        fragments.write (str(i))
        fragments.write ("\n")
        fragments.close()
        exec(open("./Search scripts/25-dir.py").read())
        fragments = open(r"25-fragments-d.txt", "a+")
        fragments.write ("\n")
        fragments.close()

for i in range (1,11): # 10 iterations
        match = open(r"30-dir.txt", "a+") #15 b matching fragments numbers
        match.write ("Iteration ")
        match.write (str(i))
        match.write ("\n")
        match.close()
        fragments = open(r"30-fragments-d.txt", "a+") #15 b matching fragments sequences and positions
        fragments.write ("Iteration ")
        fragments.write (str(i))
        fragments.write ("\n")
        fragments.close()
        exec(open("./Search scripts/30-dir.py").read())
        fragments = open(r"30-fragments-d.txt", "a+")
        fragments.write ("\n")
        fragments.close()

#For the reverse chromosome DNA chain:

for i in range (1,11): # 10 iterations
        match = open(r"5-rev.txt", "a+") #15 b matching fragments numbers
        match.write ("Iteration ")
        match.write (str(i))
        match.write ("\n")
        match.close()
        exec(open("./Search scripts/5-rev-restricted.py").read())

for i in range (1,11): # 10 iterations
        match = open(r"10-rev.txt", "a+") #15 b matching fragments numbers
        match.write ("Iteration ")
        match.write (str(i))
        match.write ("\n")
        match.close()
        exec(open("./Search scripts/10-rev-restricted.py").read())

for i in range (1,11): # 10 iterations
        match = open(r"15-rev.txt", "a+") #15 b matching fragments numbers
        match.write ("Iteration ")
        match.write (str(i))
        match.write ("\n")
        match.close()
        fragments = open(r"15-fragments-r.txt", "a+") #15 b matching fragments sequences and positions
        fragments.write ("Iteration ")
        fragments.write (str(i))
        fragments.write ("\n")
        fragments.close()
        exec(open("./Search scripts/15-rev.py").read())
        fragments = open(r"15-fragments-r.txt", "a+")
        fragments.write ("\n")
        fragments.close()

for i in range (1,11): # 10 iterations
        match = open(r"20-rev.txt", "a+") #15 b matching fragments numbers
        match.write ("Iteration ")
        match.write (str(i))
        match.write ("\n")
        match.close()
        fragments = open(r"20-fragments-r.txt", "a+") #15 b matching fragments sequences and positions
        fragments.write ("Iteration ")
        fragments.write (str(i))
        fragments.write ("\n")
        fragments.close()
        exec(open("./Search scripts/20-rev.py").read())
        fragments = open(r"20-fragments-r.txt", "a+")
        fragments.write ("\n")
        fragments.close()

for i in range (1,11): # 10 iterations
        match = open(r"25-rev.txt", "a+") #15 b matching fragments numbers
        match.write ("Iteration ")
        match.write (str(i))
        match.write ("\n")
        match.close()
        fragments = open(r"25-fragments-r.txt", "a+") #15 b matching fragments sequences and positions
        fragments.write ("Iteration ")
        fragments.write (str(i))
        fragments.write ("\n")
        fragments.close()
        exec(open("./Search scripts/25-rev.py").read())
        fragments = open(r"25-fragments-r.txt", "a+")
        fragments.write ("\n")
        fragments.close()

for i in range (1,11): # 10 iterations
        match = open(r"30-rev.txt", "a+") #15 b matching fragments numbers
        match.write ("Iteration ")
        match.write (str(i))
        match.write ("\n")
        match.close()
        fragments = open(r"30-fragments-r.txt", "a+") #15 b matching fragments sequences and positions
        fragments.write ("Iteration ")
        fragments.write (str(i))
        fragments.write ("\n")
        fragments.close()
        exec(open("./Search scripts/30-rev.py").read())
        fragments = open(r"30-fragments-r.txt", "a+")
        fragments.write ("\n")
        fragments.close()

#III. OUTPUT FILES PREPARATION FOR ANALYSIS (SORTING)
exec(open("./Sorting scripts/Sorting-initial-fragments-dir.py").read()) #Sorts all random (matching and dismatching) fragments according to their positions in segment
                                                                        #(direct DNA chain). Input file:*-dir.all.txt, output file: *-dir-all.sort.txt
exec(open("./Sorting scripts/Sorting-initial-fragments-rev.py").read()) #Does the same for reverse DNA chain
exec(open("./Sorting scripts/Sorting-matching-fragments-dir.py").read()) #Sorts matching fragments according to their positions in segment and chromosome (direct DNA chain)
                                                                         # Input file: *-d.txt, output file: *-d-list.txt *-d-count.txt
exec(open("./Sorting scripts/Sorting-matching-fragments-rev.py").read()) #Does the same for reverse DNA chain
exec(open("./Sorting scripts/Sorting-chromosome-segments-dir.py").read()) #Sorts numbers of matching fragments according to chromosome segment numbers
                                                                          #and calculates average number normalized to the total fragments number (direct DNA chain)
exec(open("./Sorting scripts/Sorting-chromosome-segments-rev.py").read()) #Does the same for reverse DNA chain

#IV. ANALYSIS
# Table with the normalized average numbers of matching segments
exec(open("./Analysis scripts/Normalization-table-dir.py").read())#Generates table of normalized fragment frequencies of each length for each chromosome segment (direct DNA strain)
exec(open("./Analysis scripts/Normalization-table-rev.py").read())#Does the same for reverse DNA chain
exec(open("./Analysis scripts/Normalization-column-dir.py").read())#Generates normalized fragment frequencies for each chromosomal disc and calculates Pearson correlation between fragment frequencies and ectopic contacts frequencies (direct DNA strain)
exec(open("./Analysis scripts/Normalization-column-rev.py").read())#Does the same for reverse DNA chain
exec(open("./Analysis scripts/Column-combinations.py").read())#Combines columns with ectopic contacts and fragment frequencies for both DNA chains
exec(open("./Analysis scripts/Normalization-column-both.py").read())#Calculates Pearson correlation between fragment frequencies and ectopic contacts frequencies (both DNA strain)
