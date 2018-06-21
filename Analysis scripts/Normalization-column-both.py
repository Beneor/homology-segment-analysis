# This script generates columns of normalized fragment frequencies for corresponding chromosomal discs
# and calculates Pearson correlation between fragment frequencies and ectopic contact frequencies
# (both DNA chains: summary fragment frequency)

# Modules to be used
import glob
import os
import math

#I. Pearson correlation calculation

fileSuffix='-column-both.txt'
filesMask='*'+fileSuffix
resultFiles = glob.glob(filesMask)

pearsonfile = open(r"Pearson-data-both.txt","a+")
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

fileSuffix='-temporal-both.txt' #Deletes -new files
filesMask='*'+fileSuffix
resultFiles = glob.glob(filesMask)
for f in resultFiles:
        os.remove(f)

