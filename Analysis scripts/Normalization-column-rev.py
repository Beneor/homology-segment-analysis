# This script generates columns of normalized fragment frequencies for corresponding chromosomal discs
# and calculates Pearson correlation between fragment frequencies and ectopic contact frequencies
# (reverse DNA chain)

# Modules to be used
import glob
import os
import math

# I. Column generation start
fileSuffix='-norm-rev.txt'
filesMask='*'+fileSuffix
resultFiles = glob.glob(filesMask)
newfileSuffix='-temporal-rev.txt'

for fileName in resultFiles:
        normalized=open(fileName,"r")
        column=open(fileName.rstrip(fileSuffix)+newfileSuffix,"w+")
        segmentnumbers = open(r"Rev-numbers.txt","r")
        for line in segmentnumbers:
                x = str(line)
                num = int(x.strip())
                bands=list()
                if num in range(1198,1278):continue #!!! SELECT THE AREA Optional: excluding segments corresponding to Segment 0 from calculations
                column.write(str(num)),column.write(";")
                normalized.seek(0)
                for line in normalized:
                        bands.append(line)
                column.write(bands[num].strip())
                column.write("\n")
        normalized.close()
        column.close()
        segmentnumbers.close()

#II.Column preparation
# 1. Frequencies assignment to discs

fileSuffix='-temporal-rev.txt'
newfileSuffix='-column-rev.txt'
filesMask='*'+fileSuffix
resultFiles = glob.glob(filesMask)

ecffile=open(r"ECF.txt","r") #Chromosome bands borders and ECF values

for fileName in resultFiles:
        column=open(fileName,"r")
        newcolumn=open(fileName.rstrip(fileSuffix)+newfileSuffix,"w+")
        newcolumn.write("Band;fragmentFrequency;ECF\n")
        ecffile.seek(0)
        for line in ecffile:
                x=str(line)
                band=x[0:x.find(",")]
                a=int(x.find(",",x.find(",")+1))
                start=int(x[x.find(",")+1:a])
                end=int(x[a+1:x.rfind(",")])
                startnum=round(start/10000)+1 #Number of the first chromosome fragment for the corresponding disc
                endnum=round(end/10000)#Number of the last chromosome fragment for the corresponding disc
                ecf=x[x.rfind(",")+1:].strip()#Ectopic contacts frequency
                newcolumn.write(band), newcolumn.write(";")
                sumfrfreq=0
                column.seek(0)
                num=0
                for line in column:
                        y=str(line)
                        frfreq=float(y[y.find(";")+1:])
                        num=int(y[0:y.find(";")])
                        if num < startnum: continue
                        if num > endnum: continue
                        sumfrfreq+=frfreq #Summary fragment frequency for corresponding disc
                        num+=1 #Summary number of bands
                dec=round(sumfrfreq,6)
                newcolumn.write("{:0.10f}".format(dec))
                newcolumn.write(";")
                newcolumn.write(ecf)
                newcolumn.write("\n")
        column.close()
        newcolumn.close()
ecffile.close()

#III. Pearson correlation calculation

fileSuffix='-column-rev.txt'
filesMask='*'+fileSuffix
resultFiles = glob.glob(filesMask)

pearsonfile = open(r"Pearson-data-rev.txt","a+")
pearsonfile.write("Filename;Pearson correlation;N;Significance\n") 
for fileName in resultFiles:
        newcolumn=open(fileName,"r")
        frfreq=0
        ecf=0
        n=0
        x1=0
        y1=0
        x4=0
        y4=0
        h2=0
        for line in newcolumn:
                z = str(line)
                if "Band" in z: continue
                frfreq = float(z[z.find(";")+1:z.rfind(";")])
                ecf = float(z[z.rfind(";")+1:])
                x1+=frfreq #Sum of fragment frequencies*ect
                y1+=ecf #Sum of fragment frequencies
                n+=1
        avfr=x1/n #average fragment frequency
        avecf=y1/n #average ectopic contact frequency
        newcolumn.seek(0)
        for line in newcolumn:
                z = str(line)
                if "Band" in z: continue
                frfreq = float(z[z.find(";")+1:z.rfind(";")])
                ecf = float(z[z.rfind(";")+1:])
                x2=frfreq-avfr
                x3=math.pow(x2,2)
                y2=ecf-avecf
                y3=math.pow(y2,2)
                h=x2*y2
                x4+=x3 #Sum of squares
                y4+=y3 #Sum of squares
                h2+=h #Sum of multiplications
        cov = h2/(math.sqrt(x4*y4))
        cov2=round(cov,3)
        pearsonfile.write(fileName),pearsonfile.write(";"),pearsonfile.write ("{:0.10f}".format(cov2)),pearsonfile.write(";")

        # Sifnificance of Pearson correlation
        tconf = cov*math.sqrt(n)/(1-math.pow(cov,2))
        #Level of confidence at alpha = 5%
        #choosing Student coefficient for n and alpha
        pearsonfile.write(str(n)),pearsonfile.write(";")
        if n >= 120:
                trep = 1.98
        elif 60 <= n <120:
                trep = 2.00
        elif 40 <=n < 120:
                trep = 2.02
        else:
                trep = 2.06
        if tconf >  trep:
                pearsonfile.write("Yes\n")
        else:
                pearsonfile.write("No\n")

# IV. Deletion of temporal files

fileSuffix='-temporal-rev.txt' #Deletes -new files
filesMask='*'+fileSuffix
resultFiles = glob.glob(filesMask)
for f in resultFiles:
        os.remove(f)

