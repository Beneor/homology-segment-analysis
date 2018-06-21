#This number sorts chromosome segments according to their number and calculates average values
# of matching fragment frequency for 10 iterations (direct DNA chain)

# Modules to be used
import glob
import os

a = "Iteration"
b = "Matching"

# Chromosomal analysis result files
fileSuffix='-dir.txt'
newfileSuffix='-sort-dir.txt'
newnewfileSuffix='-aver-dir.txt'
filesMask='*'+fileSuffix
resultFiles = glob.glob(filesMask)

# I. Sorting strings by numbers (for direct sequence)

def sortnumber(line): #Sorting lines according to the files numbers in strings
        x = str(line)
        return int(x[22:x.find("-")].strip())

for fileName in resultFiles:
        f=open(fileName, "r")  # Files with the number of matching segments - initial
        f1=open(fileName.rstrip(fileSuffix)+newfileSuffix,"w") # Sorted file
        bands = list() 
        print("%-12s"%(f.name))
        for line in f:
                if a in line: continue
                if b in line: continue
                bands.append(line) #Creating list containing all strings 
        bands.sort(key=sortnumber) #Sorting values in the list
        for i in bands:
                f1.write(i[0:i.find("-")].strip()) #Writing valuse of the list as separate strings
                f1.write("\n")
        f1.close()
        f.close()

# II. Calculating the average values for 10 iterations:

newfilesMask='*'+newfileSuffix
newresultFiles = glob.glob(newfilesMask)

for fileName in newresultFiles:
        f=open(fileName, "r")
        f2=open(fileName.rstrip(newfileSuffix)+newnewfileSuffix,"w")
        while True:
                ms = 0 #Matching segments sum
                ts = 0 #Total segments sum
                for i in range(10):
                        line = f.readline() #reading next 10 lines from input file
                        if not line:
                                break
                        line = line.strip()
                        m=float(line[0:11].strip()) #Matching segments number
                        t=float(line[12:18].strip()) #Total segments number
                        n=int(line[22:].strip()) #Chromosome part number
                        ms+=m
                        ts+=t
                        av=round(ms/ts,6) #Average number of matching segments normalized to average number of segments        
                if not line:
                        break
                f2.write("{:0.10f}  ;{:5n}\n".format(av,n))
        f2.close()
        f.close()

# III Deletion of temporal files

fileSuffix ='-sort-dir.txt' 
filesMask='*'+fileSuffix
resultFiles = glob.glob(filesMask)
for f in resultFiles:
        os.remove(f)
