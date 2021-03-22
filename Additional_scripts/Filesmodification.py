#!/usr/bin/env python3

import glob   
import os
import shutil

programDescription = '''
! Place this script into the folder 'examples' (if absent) and run it in that folder using Python3. 

This script arranges files produced by segmentanalysis in folder 'examples', placing them to folder Cytocountsdata: fragments for bands or zones: 'discs' or 'zones' folder; fragments without and with exclusion of DNA repeats: 'all' or 'nr' folder; fragments sequences or fragments frequencies  - 'l*-fragments' or 'l*-cytocounts' folder, accodring to fragments length *.
'''

a=input('Enter chromosomal areas, discs or zones: ')
b= input('Enter fragments, all or nr: ')
path='./'+a+'/'+b
os.mkdir('./'+a)
os.mkdir('./'+a+'/'+b)

s='.txt'
t='*'+s
borders={}
if a =='discs':
    f = open(r'discsfile','r')
elif a =='zones':
    f = open(r'zonesfile','r')
for line in f:
        x = str(line).strip()
        x1 = x[:x.find('\t')].strip()
        x2 = x[x.find('\t')+1:]
        borders[x1]=x2

for root, dirs, files in os.walk("Results"):
    for filename in files:
        if filename.startswith("ncounts"):
            os.remove(os.path.join(root, filename))
for root, dirs, files in os.walk("Results"):
    for filename in files:
        if filename.startswith("fragments"):            
            shutil.move (os.path.join(root, filename), (os.path.join(path, filename)))
for root, dirs, files in os.walk("Results"):
    for filename in files:
        if filename.endswith("merged.txt"):            
            shutil.move (os.path.join(root, filename), (os.path.join(path, filename)))
for root, dirs, files in os.walk("Results"):
    for filename in files:
        if filename.startswith("cytocounts"):
            os.remove(os.path.join(root, filename))
for root, dirs, files in os.walk("Results"):
    for dirname in dirs:
            os.rmdir(os.path.join(root, dirname))
os.rmdir("Results") 

cytocountdic={}
cytocountset=set()
cytocountset2=set()
fragmdic={}

prefix=''
if b=='nr':
    prefix='nr'

a= os.listdir(path)
for i in a:
    if 'director' in i:
        continue
    nametype=i[i.rfind('/')+1:i.find('.')]
    namelength=i[i.find('.')+1:i.find('-')].strip()
    name=namelength+prefix+'-'+nametype
    cytocountset.add(name)
    name2=nametype+'.'+namelength
    cytocountset2.add(name2)
for i in cytocountset:
    directory=path+'/'+i
    os.mkdir(directory)
for name2 in cytocountset2:    
    for filename in glob.glob(os.path.join(path,t)):
        if name2 in filename:
            nametype=name2[:name2.find('.')]
            namelength=name2[name2.find('.')+1:]
            name=path+'/'+namelength+'-'+nametype
            shutil.move(filename,name)


for path, dirs, files in os.walk(path):
    for filename in files:
        fullpath = os.path.join(path, filename)
        file1 = open(fullpath,'r')
        for key in borders:
            if key in filename:
                filename2 = borders[key]
                if 'fragment' in filename:
                    if 'dir' in filename:
                        filename2+='-dir'
                    if 'rev' in filename:
                        filename2+='-rev'
                filename2+=s
                fullpath2 = os.path.join(path, filename2)
                os.rename(fullpath,fullpath2)
        file1.close()
f.close()
        
    
