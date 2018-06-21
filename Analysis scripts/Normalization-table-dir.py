#This script normalizes the average fragment frequencies and produces frequency table
# (for direct DNA chain)

#!/usr/bin/env python3

# Modules to be used
import glob
import os

# Chromosomal analysis result files
fileSuffix='-aver-dir.txt'
filesMask='*'+fileSuffix
resultFiles = glob.glob(filesMask)
newfileSuffix='-norm-dir.txt'

# I. Normalization of average values

for fileName in resultFiles:
        average=open(fileName, "r")  # Files with the average number of matching segments
        normalized = open(fileName.rstrip(fileSuffix)+newfileSuffix,"w+")
        normalized.write(fileName.rstrip(fileSuffix)),normalized.write("\n")
        summa=0
        number=0
        for line in average:
                x=str(line)
                num=int(x[x.find(";")+1:].strip())
                if num in range(1198,1278):continue #!!! SELECT THE AREA! Optional: excluding segments corresponding to Segment 0 from calculations
                value=float(x[0:x.find(";")].strip())
                summa+=value
                number+=1
        average.seek(0)
        for line in average:
                x = str(line)
                value=float(x[0:x.find(";")].strip())
                norm = float(number*value/summa) #Value normalized to the average value for the chromosome
                dec = round(norm,6)
                normalized.write("{:0.10f}".format(dec))
                normalized.write("\n")
        normalized.close()



# II. Chromosome segments number output
segments=open(r"20-aver-dir.txt","r")
segmentnumbers=open(r"Dir-numbers.txt","w+")
for line in segments:
        x = str(line)
        num=int(x[x.find(";")+1:].strip())
        segmentnumbers.write(str(num)),segmentnumbers.write("\n")
segmentnumbers.close()

# III. Table generation
table = open(r"Table-dir.txt","w+")
table.write("N;10-dir;15-dir;20-dir;25-dir;30-dir;5-dir\n")
segmentnumbers = open(r"Dir-numbers.txt","r")
newfilesMask='*'+newfileSuffix
newresultFiles = glob.glob(newfilesMask)
for line in segmentnumbers:
        x = str(line)
        num = int(x.strip())
        if num in range(1198,1278):continue #!!! SELECT THE AREA Optional: excluding segments corresponding to Segment 0 from calculations
        table.write(str(round(num,0))),table.write(";")
        for fileName in newresultFiles:
                temporal=open(fileName,"r")
                bands=list()
                for line in temporal:
                        bands.append(line)
                table.write(bands[num].strip()),table.write(";")
                temporal.close()
        table.write("\n")

