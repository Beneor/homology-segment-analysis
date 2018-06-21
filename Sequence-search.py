import glob
import math

#I. Segment
fileSuffix='-fr.txt'
filesMask='*'+fileSuffix
resultFiles = glob.glob(filesMask)
seqsearch=open(r"Sequence-search.txt","w+")
x="aattttcgcattttttgtaaggggtaacatcatcaaaatttgcaaaaaat"# INPUT THE SEQUENCE!!!

# 2. Preparation of complementary segment
c=""
for i in range(0,len(x)):
    j = int(len(x)-i-1)
    s = x[j]
    if s == "a":
            c+="t"
    elif s == "g":
            c+="c"
    elif s == "c":
            c+="g"
    elif s == "t":
            c+="a"
    else:
            c+="n"

#III. Segments search
for fileName in resultFiles:
    f=open(fileName, "r")  #Chromosomal parts files
    for line in f:
            a=str(line)
            p=a.count(x) #Number of segment in direct DNA strain
            q=a.count(c) #Number of segment in reverse DNA strain
            summary=p+q #Summary number in both strains
    print("%-12s%-17s"%(summary,f.name))
    seqsearch.write(str(summary))
    seqsearch.write(";")
    seqsearch.write(f.name)
    seqsearch.write("\n")
    f.close()
seqsearch.close()

#II. Sorting chromosome segments
def sortnumber(line): #Sorting lines according to the files numbers in strings
        e = str(line)
        return int(e[e.rfind(";")+1:e.find("-")].strip())#The first segment position
seqsearch=open(r"Sequence-search.txt","r")  # Initial order
sortedMatchingFragments=open(r"Sequence-sort.txt","w+") # Sorted order
bands=list()
for line in seqsearch:
    bands.append(line) #Creating list containing all strings
    bands.sort(key=sortnumber) #Sorting values in the list
for i in bands:
    sortedMatchingFragments.write(i) #Writing valuse of the list as separate strings
sortedMatchingFragments.close()
seqsearch.close()

#III.Column preparation
# 1. Frequencies assignment to discs
ecffile=open(r"ECF.txt","r") #Chromosome bands borders and ECF values
column=open(r"Sequence-sort.txt","r")
newcolumn=open(r"Sequence-sort-ecf.txt","w+")
newcolumn.write("Band;fragmentFrequency;ECF\n")
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
            frfreq=float(y[0:y.find(";")])
            num=int(y[y.find(";")+1:y.find("-")])
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

#IV. Pearson correlation calculation

newcolumn=open(r"Sequence-sort-ecf.txt","r")
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
print("\n")
print("Pearson correlation with ECF:")
print("Covariation coefficient: "), print(cov) #Pearson correlation

# Reliability of Pearson correlation

tconf = cov*math.sqrt(n)/(1-math.pow(cov,2))
print("t statistics: ")
print(tconf)
#Level of confidence at alpha = 5%
#choosing Student coefficient for n and alpha
if n >= 120:
    trep = 1.98
elif 60 <= n <120:
    trep = 2.00
elif 40 <=n < 120:
    trep = 2.02
else:
    trep = 2.06
if tconf >  trep: 
    print("Reliable at p<0.05")
else:
    print("Unreliable at p<0.05")

