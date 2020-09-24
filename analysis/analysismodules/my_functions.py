#!/usr/bin/env python3

import glob
import os
from scipy.stats.stats import spearmanr
from scipy import stats
import shutil
import math

# Accessory fuctions used for analysis of Homology Segment Analysis results

def averagecorr(fragname,length,s):
    ''' Returns discs with significant correlations, calculates the average significant Spearman's rank correlation and the share of significant probabilities'''
    result = open(fragname+'-'+str(length)+'-analysis-average-correlation.xls','w') #Result file vith average data and standard errors
    filesuffix=s
    filemask='*'+filesuffix
    files=glob.glob(filemask)
    corvoc={} #Significant correlations and probabilities data for the given fragment length
    discvoc={} #Discs with significant correlations for the given fragment length
    for f1 in files: #Files of correlations (csv2) for the given fragment lengths
        f=open(f1,'r')
        v=0 # p <P
        e=0 # no ectopics
        n=0 # total number
        j=0 #summary correlation
        length = int(f1[:f1.find('-')]) #Name of the file (length)
        corsigniflist=[] #Significant correlations
        probsigniflist=[] #Significant correlation probabilities
        discsignlist=[]#Discs with significant correlations
        for line in f:
            n+=1
            x=str(line)
            if 'nan' in x:
                e+=1
            else:
                signif = float(x[x.rfind("\t")+1:].strip()) #Correlation significance level
                if x.count('\t')==2:
                    cor = float(x[x.find("\t")+1:x.rfind("\t")])
                    dsk=x[:x.find('\t')]
                elif x.count('\t')==3:
                    cor0= x[x.find('\t')+1:]
                    cor = float(cor0[cor0.find('\t')+1:cor0.rfind('\t')])
                    dsk0=x[:x.rfind('\t')]
                    dsk=dsk0[:dsk0.rfind('\t')]
                else:
                    print ("Non-standard correlation file\n")
                if signif <0.05:
                    corsigniflist.append(cor)
                    probsigniflist.append(signif)
                    v+=1 #Number of significant correlations
                    j+=cor
                    discsignlist.append(dsk)
        corsterr=str(round(stats.sem(corsigniflist),3))  #Correlations standard error calculations
        probsterr=str(round(stats.sem(probsigniflist),3))  #Probabilities standard error calculations     
        h=n-e #Total number of correlations
        if h!=0:
            a1= v/h # The percentage of significant correlations
        else:
            a1=0
        if v!=0:
            avcor = str(round((j/v),3))
        else:
            avcor='no'
        dat=avcor+'\t'+str(round(a1,3))+'\t'+str(h)+'\t'+str(v)+'\t'+corsterr+'\t'+probsterr
        corvoc[length]=dat
        discvoc[length]=discsignlist
        f.close()
    listcor = list(corvoc.keys())
    listcor.sort()
    listdiscs = list(discvoc.keys())
    listdiscs.sort()
    result.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%('length','rho', 'p','total number of correlations','number of significant correlations','rho_sterr','p_sterr'))
    for i in listcor:
        result.write('%s\t%s\n'%(i,corvoc[i]))
    result.close()
    for i in listdiscs:
        resultdiscs=open(str(i)+'l-analysis-discs','w') #File containing discs with significant correlations
        for j in discvoc[i]:
            resultdiscs.write('%s\t%s\n'%(i,j))
    resultdiscs.close()
    
def combs(n):
    """Returns list of lists of all combinations ofj in i 1 and 0 for n values"""
    import copy
    a=[[]]
    b=copy.deepcopy(a)
    for k in range(n):
        for i in range(len(a)):
            a[i].append(0)
        for j in range(len(b)):
            b[j].append(1)
        c=a+b
        a=c
        b=copy.deepcopy(a)
    return(c)

def comparetup (t1,t2):
    """Compares two tuples and returns the number of mismatches"""
    n=0
    for i in range(len(t1)):
       a=t1[i]
       b=t2[i]
       if a !=b:
           n+=1
    return(n)

def disccount(length):
    """Returns disc pairs with significant unspecific Spearman correlation"""
    a=[]
    b=[]
    file1=open(str(length)+'l-analysis-discs','r')
    file2=open(str(length)+'l-discpairs-unspecific-counts.csv3','w')
    for line in file1:
        x = str(line)
        name = int(x[:x.find('\t')])
        frfd=int(x[x.find('\t')+1:x.rfind('\t')])
        ecfd=int(x[x.rfind('\t')+1:])
        a.append(frfd)
        b.append(ecfd)
    c = set(a)
    d = set(b)
    e=list(c)
    e.sort()
    f=list(d)
    f.sort()
    
    def discorarr(p,q):
        """Arrange discs according to number of their correlations"""
        voc={}
        for i in q:
            j = p.count(i)
            voc[i] = j 
        l=list(voc.keys())
        l.sort() #Keys - the sorted values of disc numbers
        for k in l:
            file2.write('%s\t%s\n'%(k, voc[k])) #First - disc number, second - disc frequency
    file2.write('FRF discs:\n')
    discorarr(a,e)
    file2.write('ECF discs:\n')
    discorarr(b,f)
    
def discequal(length):
    """Returns discs with significan specific Spearman correlation"""
    file1=open(str(length)+'-specifcor.csv2','r')
    file2=open(str(length)+'l-disc-specific-counts.csv3','w')
    for line in file1:
        x = str(line)
        name = int(x[:x.find('\t')])
        if 'nan' in x:
            file2.write('%s\t%s\n'%(name,'nan'))
            continue
        corr = float(x[x.find('\t')+1:x.rfind('\t')])
        prob=float(x[x.rfind('\t')+1:].strip())
        if prob >= 0.05:
            file2.write('%s\n'%(name))
        if prob < 0.05:
            file2.write('%s\t%f\t%f\n'%(name,corr,prob))
    file1.close()
    file2.close()
        
def mannwhitney(path1,path2,u,name1,name2,x,ex):
    """Calculates Mann-Whitney U statistics for Spearman rho and p files from path1 = p1 and path2 = p2. Also tests samples for normality usinf Shapiro Wilk test"""
    result=open (x+'-'+name1+'-'+name2+'-MW-test.xls','w')
    result.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%('name1','name2','prob_cor','prob_p','cordistr1','cordistr2','probdistr1','probdistr2'))
    def voc(a,b):
        """ Returns dictionary containing filenames as keys and lists of discs (c,d) with Spearman correlation rho coefficients and probabilities. Also checks the arranged list of keys (e)"""
        c={} #Dictionary for filenames (keys) and lists of rho values)
        d={} #Dictionary for filenames (keys) and lists of p values)
        e=[] #List for the arranged keys in c and d 
        for filename in glob.glob(os.path.join(a,b)):
            f1 = open(filename,'r')
            fn=filename[filename.rfind('/')+1:filename.rfind('-')] #File number
            a=[]
            b=[]
            for line in f1:
                x=str(line)
                if 'nan' in x:
                    continue
                prob=float(x[x.rfind('\t')+1:].strip())
                if x.count('\t')==2:
                    cor=float(x[x.find('\t')+1:x.rfind('\t')])
                elif x.count('\t')==3:
                    cor0=x[:x.rfind('\t')]
                    cor=float(cor0[cor0.rfind('\t')+1:].strip())
                a.append(cor)
                b.append(prob)
            c[fn]=a #correlation value
            d[fn]=b #correlation probability
            e.append(fn)
            f1.close()
        e.sort()
        return c, d, e
    
    c1,c2,e = voc(a=path1,b=ex)
    d1,d2,g = voc(a=path2,b=ex)

    for k in range(len(e)):
        w1=[]
        i=e[k]
        for m in range(len(g)):
            j=g[m]
            if u == True: #Parameter that specifies whether files are from the same samplings
                if i==j:
                    continue
            else:
                if i!=j:
                    continue
            cor1=c1[i] 
            cor2=d1[j]
            prob1=c2[i]
            prob2=d2[j]
            corv,corp = stats.mannwhitneyu(cor1,cor2,alternative='two-sided') #Two sided Mann-Whitney U test
            c1v,c1p =stats.shapiro(cor1)
            c2v,c2p =stats.shapiro(cor2)
            p1v,p1p =stats.shapiro(prob1)
            p2v,p2p =stats.shapiro(prob2)
            cordist1, cordist2, probdist1, probdist2 = 'norm' #Variables for samples distributions
            if c1p <0.05:
                cordist1='not-norm'
            else:
                cordist1='norm'
            if c2p <0.05:
                cordist2='not-norm'
            else:
                cordist2='norm'
            if p1p <0.05:
                probdist1='not-norm'
            else:
                probdist='norm'
            if p2p <0.05:
                probdist2='not-norm'
            else:
                probdist2='norm'
            probv,probp = stats.mannwhitneyu(prob1,prob2,alternative='two-sided')
            corp1=round(corp,3) #Mann-Whitney U test; p value for rho correlations comparision 
            probp1=round(probp,3) #Mann-Whitney U test; p value for rho probabilities comparision
            if corp1 <0.05 or probp1 <0.05:
                result.write('%s\t%s\t%f\t%f\t%s\t%s\t%s\t%s\n'%(i,j,corp1,probp1,cordist1,cordist2,probdist1,probdist2))
    result.close()

def movefiles (a):
    ''' Copies files to the working directory'''
    for root, dirs, files in os.walk(a):
        for file in files:
            shutil.copy(os.path.join(root, file), (os.path.join(r'../Analysis', file)))
            
def nearest(n,length,fragdir,ectdir):
    """ Sequentially calculates Spearman's rank correlation for n nearest zones"""
    p = '../Analysis' #Path to FRF files
    result = open(str(n)+'-'+str(length)+'-nearest.csv2','w')
    def freqdic(path,b,t):
        """ Returns the dictionary (c) with keys for the filenames and values for massive of frequencies for each line within a file. Also returns the sorted list (l) for the integer representations of filenames"""
        g=[[]]
        c={}
        i=0
        e = "X\t"
        for filename in glob.glob(os.path.join(path,b)):
            f = open(filename,'r')
            f1=int(filename[filename.rfind('/')+1:].strip(b))
            for line in f:
                x=str(line)
                if t is True and e not in x:
                    continue
                x1 =float(x[x.rfind('\t')+1:].strip())
                g[i].append(x1)
            g.append([])
            c[f1]=g[i]
            f.close()
            i+=1
        l=list(c.keys())
        l.sort()
        return c, l

    # I. Generation of dictionaries and lists for FRF and ECF
    c1,l1 = freqdic(path=fragdir,b='*.txt',t=True)
    c2,l2 = freqdic(path=ectdir,b='*.csv',t=False)

    # II. Generation of FRF and ECF samples for n file pairs containing N FRF and ECF data. # Spearman rank coefficient correlation for each FFR and ECF samples and data output
    for k in range (len(l1)-n+1):
        w1=[] # FRF samaple
        w2=[] # ECF sample
        i = l1[k] # The name of FRF file, arranged by the order
        j = l2[k] # The name of ECF file, arranged by the order
        if i != j:   # FRF and ECF files numbers should be equal (control)
            continue
        for z in range(k,k+n):
            x=float(c1[i][z]) # FRF value for the z-th disc (z is between k and k+n).
            y=float(c2[j][z]) # FRF value for the i-th disc in the i-th file
            w1.append(x) #FRF sample
            w2.append(y) #ECF sample
        rho, pval = spearmanr(w1,w2)
        result.write('%s\t%f\t%f\n'%(l1[k],round(rho,3), round(pval,3)))
    result.close()

def pathgenerator(a,t): 
    """Generates path from the list of directories separated by t"""
    j=''
    for i in a:
        j+=(i+t)
    return(j.rstrip(t))

def rename(t,s):
    '''Renames txt files according to list'''
    filesuffix = s 
    filemask = '*'+ filesuffix
    files = glob.glob(filemask)

    a={}
    discs = open(t,'r') # Opens the list containing names of files
    for line in discs:
        x = str(line).strip()
        x1 = x[:x.find('\t')].strip()
        x2 = x[x.find('\t')+1:]
        a[x1]=x2
    for filename in files:
        file1 = open(filename,'r')
        for key in a:
            if key in filename:
                filename2 = a[key]+s
                file2 = open(filename2,'w')
                for line in file1:
                    x = str(line)
                    file2.write(x)
                file2.close()
        file1.close()
        os.remove(filename)
    discs.close()

def removefiles(s):
    '''Removes files with specified ending (s) in the working directory'''
    filesuffix = s 
    filemask = '*'+ filesuffix
    files = glob.glob(filemask)
    for filename in files:
        os.remove(filename)    
        
def specifcor(length,fragdir,ectdir):
    """ Sequentially calculates Spearman's rank correlation for n nearest zones"""
    p = '../Analysis' #Path to FRF files
    result = open(str(length)+'-specifcor.csv2','w')
    def freqdic(path,b,t):
        """ Returns the dictionary (c) with keys for the filenames and values for massive of frequencies for each line within a file. Also returns the sorted list (l) for the integer representations of filenames"""
        g=[[]]
        c={}
        i=0
        e = "X\t"
        for filename in glob.glob(os.path.join(path,b)):
            f = open(filename,'r')
            f1=int(filename[filename.rfind('/')+1:].strip(b))
            for line in f:
                x=str(line)
                if t is True and e not in x:
                    continue
                x1 =float(x[x.rfind('\t')+1:].strip())
                g[i].append(x1)
            g.append([])
            c[f1]=g[i]
            f.close()
            i+=1
        l=list(c.keys())
        l.sort()
        return c, l

    # I. Generation of dictionaries and lists for FRF and ECF
    c1,l1 = freqdic(path=fragdir,b='*.txt',t=True)
    c2,l2 = freqdic(path=ectdir,b='*.csv',t=False)

    # II. Generation of FRF and ECF samples for n file pairs containing N FRF and ECF data. # Spearman rank coefficient correlation for each FFR and ECF samples and data output
    for k in range (len(l1)):
        w1=[] # FRF samaple
        w2=[] # ECF sample
        i = l1[k] # The name of FRF file, arranged by the order
        j = l2[k] # The name of ECF file, arranged by the order
        if i != j:   # FRF and ECF files numbers should be equal (control)
            continue
        for z in range(len(l1)+1):
            x=float(c1[i][z]) # FRF value for the z-th disc (z is between k and k+n).
            y=float(c2[j][z]) # FRF value for the i-th disc in the i-th file
            w1.append(x) #FRF sample
            w2.append(y) #ECF sample
        rho, pval = spearmanr(w1,w2)
        result.write('%s\t%f\t%f\n'%(l1[k],round(rho,3), round(pval,3)))
    result.close()

def unspecifcor(length,fragdir,ectdir):
    """ Sequentially calculates Spearman's rank correlation for n nearest zones"""
    p = '../Analysis' #Path to FRF files
    result = open(str(length)+'-unspecifcor.csv2','w')
    def freqdic(path,b,t):
        """ Returns the dictionary (c) with keys for the filenames and values for massive of frequencies for each line within a file. Also returns the sorted list (l) for the integer representations of filenames"""
        g=[[]]
        c={}
        i=0
        e = "X\t"
        for filename in glob.glob(os.path.join(path,b)):
            f = open(filename,'r')
            f1=int(filename[filename.rfind('/')+1:].strip(b))
            for line in f:
                x=str(line)
                if t is True and e not in x:
                    continue
                x1 =float(x[x.rfind('\t')+1:].strip())
                g[i].append(x1)
            g.append([])
            c[f1]=g[i]
            f.close()
            i+=1
        l=list(c.keys())
        l.sort()
        return c, l
    # I. Generation of dictionaries and lists for FRF and ECF
    c1,l1 = freqdic(path=fragdir,b='*.txt',t=True)
    c2,l2 = freqdic(path=ectdir,b='*.csv',t=False)
    
    # II. Generation of FRF and ECF samples for file pairs containing FRF and ECF data 
    # (for different discs). 
    # Spearman rank coefficient correlation for each FFR and ECF samples and data output
    for k in range (len(l1)):
        w1=[] # FRF sample
        i = l1[k] # The name of FRF file, arranged by the order
        for z in range(len(l1)+1):
            x=float(c1[i][z]) # FRF value for the z-th disc (z is between k and k+n).
            w1.append(x) #FRF sample
        for m in range (len(l2)):
            w2=[] # ECF sample
            j = l2[m] # The name of ECF file, arranged by the order
            if i == j:   # FRF and ECF files numbers should be equal (control)
                continue
            for z in range(len(l1)+1):
                y=float(c2[j][z]) # FRF value for the i-th disc in the i-th file
                w2.append(y) #ECF sample
            rho, pval = spearmanr(w1,w2)
            result.write('%s\t%s\t%f\t%f\n'%(l1[k],l2[m], round(rho,3), round(pval,3)))
    result.close()
