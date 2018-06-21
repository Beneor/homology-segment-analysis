# This script combines the normilized fragment frequencies
# for direct and reverse DNA chains


# Modules to be used
import glob
import os

# Chromosomal analysis result files
fileSuffix='-column-dir.txt'
newfileSuffix='-column-rev.txt'
newfileSuffix2='-column-both.txt'
filesMask='*'+fileSuffix
resultFiles = glob.glob(filesMask)

a = "Band"

# I. Normalization of average values

for fileName in resultFiles:
        direct=open(fileName, "r")  # Column file for direct DNA chain
        reverse=open(fileName.rstrip(fileSuffix)+newfileSuffix,"r") # Column file for reverse DNA chain
        combined=open(fileName.rstrip(fileSuffix)+newfileSuffix2,"w+") # Column file for both DNA chains
        combined.write("Band;fragmentFrequency;ECF\n")
        band1=list()
        for line in direct:
                x = str(line)
                if a in x:
                        continue
                banddir=x[0:x.find(";")]
                freqdir=float(x[x.find(";")+1:x.rfind(";")])
                ecf=x[x.rfind(";")+1:]
                for line in reverse:
                        y = str(line)
                        if a in y: continue
                        bandrev=y[0:y.find(";")]
                        freqrev=float(y[y.find(";")+1:y.rfind(";")])
                        if bandrev == banddir:
                                comb=freqdir+freqrev
                                combined.write(banddir)
                                combined.write(";")
                                combined.write(str(comb))
                                combined.write(";")
                                combined.write(ecf)
                                break
                                
                              
                              
                              
                              
                              
                              
                              

   
