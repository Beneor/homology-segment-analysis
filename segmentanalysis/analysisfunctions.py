#!/usr/bin/env python3

# Fuctions used for analysis of Homology Segment Analysis results

import segmentanalysis.my_functions as m
import os, shutil
import glob
import copy

def dirsgen(fragname,x):
    """Generates directories for files of analysis and moves files to specific directories"""
    subd0=x+fragname
    subd1=x+fragname+'/Archive'
    subd2=x+fragname+'/Discs'
    os.mkdir(subd0)
    os.mkdir(subd1)
    os.mkdir(subd2)
    
def movefiles(fragname,x):
    """Moves files of analysis to specific directories"""
    for f in glob.glob('*'+'.csv2'): #Moves files of analysis to specific directories
        shutil.move(f,x+fragname+'/Archive') 
    for f in glob.glob('*'+'l-analysis-discs'):
        shutil.move(f,x+fragname+'/Discs')
    for f in glob.glob('*'+'.xls'):
        shutil.move(f,x+fragname)
    for f in glob.glob('*'+'.csv3'):
        shutil.move(f,x+fragname+'/Discs')
        
def correlation(strain,areas,fragments,corr,lengths,fragmentsdir,ectopicsdir,x='./Correlations/'):
    """Calculates specific and unspecific correlation for fragment and ectopic frequencies"""
    namelist=[strain,areas,fragments,corr] #Generation folder name
    fragname=m.pathgenerator(namelist,'-')
    dirsgen(fragname,x)
    for length in lengths:
        fragmlist=[fragmentsdir,areas,fragments] # Generation of the path to fragments frequencies files:
        fragdir=m.pathgenerator(fragmlist,'/')
        if fragments=='all': 
            fragdir+=('/l'+str(length)+'-cytocounts')
        elif fragments=='nr': 
            fragdir+=('/l'+str(length)+'nr-cytocounts')
        ectlist=[ectopicsdir,strain,areas] # Generation of the path to ectopic frequencies files:
        ectdir=m.pathgenerator(ectlist,'/')
        
        if corr =='specific': #Specific correlation
            m.specifcor(length,fragdir,ectdir) 
            m.discequal(length) #Discs with specific correlations
        elif corr =='unspecific': #Unspecific correlation
            m.unspecifcor(length,fragdir,ectdir)
    m.averagecorr(fragname,length,s='.csv2') #For significant correlations: calculates average Spearman correlation (rho), proportion of significant correlations (p), rho and p standard errors
    if corr=='unspecific':
        for length in lengths:
            m.disccount(length) #Disc pairs with significant unspecific correlation arranged according to numbers of correlations
    movefiles(fragname,x)
        
def neardisccor(strain,areas,fragments,corr,lengths,fragmentsdir,ectopicsdir,x='./Nearest/'):
    """Calculates specific correlation for fragment and ectopic frequencies"""
    namelist=[strain,areas,fragments,corr] #Generation folder name
    fragname=m.pathgenerator(namelist,'-')
    dirsgen(fragname,x)
    for length in lengths:
        fragmlist=[fragmentsdir,areas,fragments] # Generation of the path to fragments frequencies files:
        fragdir=m.pathgenerator(fragmlist,'/')
        if fragments=='all': 
            fragdir+=('/l'+str(length)+'-cytocounts')
        elif fragments=='nr': 
            fragdir+=('/l'+str(length)+'nr-cytocounts')
        ectlist=[ectopicsdir,strain,areas] # Generation of the path to ectopic frequencies files:
        ectdir=m.pathgenerator(ectlist,'/')
        result = open(str(length)+'-nearest-av-specific.xls','w')
        result.write('%s\t%s\t%s\n'%('area','average_correlation','average_share'))
        for n in range(10,117):
            avrholist, avproblist = m.nearest(n,length,fragdir,ectdir)
            if len(avrholist) !=0:
                avrho = sum(avrholist)/float(len(avrholist))#The average rho value for len(l1)-n chromosomal parts and n discs
            else:
                avrho = 0
            if len(avproblist) !=0:
                avpnum = sum(avproblist)/float(len(avproblist))#THe average rho value for len(l1)-n chromosomal parts and n discs
            else:
                avpnum = 0
            result.write('%s\t%f\t%f\n'%(n,round(avrho,3), round(avpnum,3)))
        result.close()
        m.removefiles('l-analysis-discs')
        subd = x+fragname+'/Archive/'+str(length) 
        os.mkdir(subd)
        for f in glob.glob('*'+'.csv2'):
            shutil.move(f,subd)
        movefiles(fragname,x)
    
def neardiscunspecifcor(strain,areas,fragments,corr,lengths,fragmentsdir,ectopicsdir,x='./Nearest/'):
    """Calculates unspecific correlation for fragment and ectopic frequencies"""
    namelist=[strain,areas,fragments,corr] #Generation folder name
    fragname=m.pathgenerator(namelist,'-')
    dirsgen(fragname,x)
    for length in lengths:
        fragmlist=[fragmentsdir,areas,fragments] # Generation of the path to fragments frequencies files:
        fragdir=m.pathgenerator(fragmlist,'/')
        if fragments=='all': 
            fragdir+=('/l'+str(length)+'-cytocounts')
        elif fragments=='nr': 
            fragdir+=('/l'+str(length)+'nr-cytocounts')
        ectlist=[ectopicsdir,strain,areas] # Generation of the path to ectopic frequencies files:
        ectdir=m.pathgenerator(ectlist,'/')
        result = open(str(length)+'-nearest-av-unspecific.xls','w')
        result.write('%s\t%s\t%s\n'%('area','average_correlation','average_share'))
        for n in range(10,117):
            avrholist, avproblist = m.nearestunspecifcor(n,length,fragdir,ectdir)
            if len(avrholist) !=0:
                avrho = sum(avrholist)/float(len(avrholist))#THe average rho value for len(l1)-n chromosomal parts and n discs
            else:
                avrho = 0
            if len(avproblist) !=0:
                avpnum = sum(avproblist)/float(len(avproblist))#THe average rho value for len(l1)-n chromosomal parts and n discs
            else:
                avpnum = 0
            result.write('%s\t%f\t%f\n'%(n,round(avrho,3), round(avpnum,3)))
        result.close()
        m.removefiles('l-analysis-discs')
        subd = x+fragname+'/Archive/'+str(length) 
        os.mkdir(subd)
        for f in glob.glob('*'+'.csv2'):
            shutil.move(f,subd)
        movefiles(fragname,x)
    
def nearestcompare(strain1,areas1,fragments1,corr1,strain2,areas2,fragments2,corr2,n,x): 
    """Mann Whitney U statistics calculation"""
    namelist1=[strain1,areas1,fragments1,corr1]
    name1=m.pathgenerator(namelist1,'-') #The name of file1
    namelist2=[strain2,areas2,fragments2,corr2]
    name2=m.pathgenerator(namelist2,'-') #The name of file1
    
    for root, dirs, files in os.walk(x):
        for dir_name in dirs:
            if name1 in dir_name:
                path1 = './'+x+'/'+name1+'/Archive/'+str(n) #x - the main directory containing files for analysis; name1 - the name of the directiry containing files set 1
    for root, dirs, files in os.walk(x):
        for dir_name in dirs:
            if name2 in dir_name:
                path2 = './'+x+'/'+name2+'/Archive/'+str(n) #name2 - the name of the directiry containing files set 1
                os.listdir(path2)
    m.nearestavcomp(path1,path2,name1,name2,n,x,ex='*-nearest.csv2') #u is True if name1==name2

def parcomb(n,dic):
    """Generates combinations of n parameters for analysis"""
    comblist=m.combs(n) #Generation of all combinations of n '0' or '1'
    for i in comblist:
        parlist=[]
        for j in range(len(i)):
            key = str(i[j])
            pardict=dic[j]
            a=pardict[key]
            parlist.append(a)
        t=tuple(parlist) #Generation of all combinations of parameters
        yield(t) #Returns combinations of parameters one by one

def parcombmdiff(n,dic,v=1):
    """Returns the tuple of unique combinations of two tuples containing n parameters with exactly m difference between two tuples in a pair """
    a=[]
    b=[]
    dir0=[]
    dir1=[]
    for t1 in parcomb(n,dic):
        a.append(t1)
    for t2 in parcomb(n,dic):
        b.append(t2)
    for i in range(len(a)):
        x=a[i]
        for j in range(len(b)):
            dir0=[]
            y=b[j]
            dir0.append(x)
            dir0.append(y)
            dir1.append(dir0)
    dir2=[]
    for k in dir1:
        x = m.comparetup(k[0],k[1])
        if x ==v:
            dir2.append(k)
    dir1=copy.deepcopy(dir2)
    for k in dir1:  #Excluding the repeated combinations for t1 and t2
        for l in dir2:
            if l[0]==k[1] and k[0]==l[1]:
                dir1.pop(dir2.index(l))
                dir2.remove(l)
    a=[]
    
    for k in dir2:
        x=k[0]+k[1]
        t=tuple(x)
        yield(t)
        
def probcompare(strain1,areas1,fragments1,corr1,strain2,areas2,fragments2,corr2,x,u):
    """Mann Whitney U statistics calculation"""
    namelist1=[strain1,areas1,fragments1,corr1]
    name1=m.pathgenerator(namelist1,'-') #The name of file1
    namelist2=[strain2,areas2,fragments2,corr2]
    name2=m.pathgenerator(namelist2,'-') #The name of file1
    
    for root, dirs, files in os.walk(x):
        for dir_name in dirs:
            if name1 in dir_name:
                path1 = './'+x+'/'+name1+'/Archive' #x - the main directory containing files for analysis; name1 - the name of the directiry containing files set 1
    for root, dirs, files in os.walk(x):
        for dir_name in dirs:
            if name2 in dir_name:
                path2 = './'+x+'/'+name2+'/Archive' #name2 - the name of the directiry containing files set 1
                os.listdir(path2)
    m.probcomp(path1,path2,name1,name2,x,ex='*.csv2') #u is True if name1==name2

    
def statmw(strain1,areas1,fragments1,corr1,strain2,areas2,fragments2,corr2,x,u): 
    """Mann Whitney U statistics calculation"""
    namelist1=[strain1,areas1,fragments1,corr1]
    name1=m.pathgenerator(namelist1,'-') #The name of file1
    namelist2=[strain2,areas2,fragments2,corr2]
    name2=m.pathgenerator(namelist2,'-') #The name of file1
    
    for root, dirs, files in os.walk(x):
        for dir_name in dirs:
            if name1 in dir_name:
                path1 = './'+x+'/'+name1+'/Archive' #x - the main directory containing files for analysis; name1 - the name of the directiry containing files set 1
    for root, dirs, files in os.walk(x):
        for dir_name in dirs:
            if name2 in dir_name:
                path2 = './'+x+'/'+name2+'/Archive' #name2 - the name of the directiry containing files set 1
                os.listdir(path2)
    m.mannwhitney(path1,path2,u,name1,name2,x,ex='*.csv2') #u is True if name1==name2

