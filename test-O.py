#!/usr/bin/python

import re
import numpy as np
import random

capList=['Neu','Fuc']
addList=['Gala3','Galb3','Galb4','GlcNAcb6','GlcNAcb3','GlcNAca6']
#addList=['Gal','Gal','GlcNAc','GlcNAc','GlcNAc']

#f=open('probabilities.dat')
f=open('conP.dat')
lines=f.readlines()
f.close()

pairDict={}
for l in lines:
    cols=l.split()
    pairDict[cols[0]]=float(cols[1])

core=np.random.choice(['Galb3GalNAc','Galb3(GlcNAcb6)GalNAc','GlcNAcb3GalNAc','GlcNAcb3(GlcNAcb6)GalNAc','GalNAca3GalNAc','GalNAcb6GalNAc','GalNAca6GalNAc','Gala3GalNAc'],p=[0.2,0.2,0.1,0.1,0.1,0.1,0.1,0.1])
base='GalNAc'

newGlycans=[]

branches=re.findall('\(.*?\)',core)
branches.append(re.sub('\(.*?\)','',core))
if core=='Galb3GalNAc':
    branches=[core]
print(core)

for branch in branches:
    branch=re.sub('[\)|\(]','',branch)
    print(branch)
    tmp=branch
    nrRes=re.split('[a|b]\d',branch)[1]
    rRes=re.split('[a|b]\d',branch)[0]
    bond0=re.findall('[a|b]\d',branch)[0]

    counter=0
    broken=False
    while counter<=3:
        nrRes=rRes
        for pair in pairDict.keys():
            if re.split('[\?|a|b][\d|\?]',pair)[1]==nrRes:
                p=pairDict[pair]
                if p>=random.random():
                    rRes=re.split('[\?|a|b][\d|\?]',pair)[0]
                    bond0=re.findall('[\?|a|b][\d|\?]',pair)[0]
                    if '%s%s'%(rRes,bond0) not in addList and rRes not in ['Fuc','Neu']:
                    #if rRes not in addList and rRes not in capList:
                        next
                    else:
                        tmp=('%s%s'%(rRes,bond0))+tmp
                        counter+=1
                        if rRes in capList:
                            broken=True
                            break
        if broken:
            break
    newGlycans.append(tmp)

if len(newGlycans)>1:
    glycan=re.sub(base,'',tmp)+'['+newGlycans[0]+']'+base
else:
    glycan=newGlycans[0]
glycan=re.sub('Neu','Neu5Ac',glycan)
a=glycan
numbers=re.compile(r'([a|b]\d)')
a=numbers.sub(r'(\1)',a)
print(a)

