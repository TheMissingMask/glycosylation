#!/usr/bin/python

import re
import random

def checkN(pair):
    rRes=re.split('[\?|a|b][\d|\?]',pair)[0]
    if rRes not in ['Man','Gal','GlcNAc','GalNAc','Fuc','Neu']:
        return False
    else:
        return True

capList=['Neu','a2','a3','a6','a4','Fuc']
addList=['Gal','GalNAc','GlcNAc','Fuc','Neu']

f=open('probabilities.dat')
lines=f.readlines()
f.close()

pairDict={}
for l in lines:
    cols=l.split()
    pairDict[cols[0]]=float(cols[1])

core='GlcNAcb2Mana6(GlcNAcb2Mana3)Manb4GlcNAcb4GlcNAc'
base=core.split(')')[1]

newGlycans=[]

branches=re.findall('\(.*?\)',core)
branches.append(re.sub('\(.*?\)','',core))


for branch in branches:

    branch=re.sub('[\)|\(]','',branch)
    tmp=branch
    nrRes=re.split('[a|b]\d',branch)[1]
    rRes=re.split('[a|b]\d',branch)[0]
    bond0=re.findall('[a|b]\d',branch)[0]

    counter=0
    broken=False
    while counter<5:
        for pair in pairDict.keys():
            if checkN(pair)==False:
                next
            if re.split('[\?|a|b][\d|\?]',pair)[1]==nrRes:
                p=pairDict[pair]
                if p>=random.random():
                    rRes=re.split('[\?|a|b][\d|\?]',pair)[0]
                    bond0=re.findall('[\?|a|b][\d|\?]',pair)[0]
                    if '?' in bond0:
                        next
                    elif rRes not in addList:
                        next
                    else:
                        tmp=('%s%s'%(rRes,bond0))+tmp
                        counter+=1
                        nrRes=rRes
                        if rRes in capList or bond0 in capList:
                            broken=True
                            break
        if broken:
            break
    newGlycans.append(tmp)

cpx=False
hybrid=False
man=False

if cpx:
    glycan=re.sub(base,'',tmp)+'['+newGlycans[0]+']'+base
elif hybrid:
    glycan='Mana6[Mana3]Mana6[%s]'%(newGlycans[0])+base
elif man:
    glycan=random.choice([
        'Mana6[Mana3]Mana6[Mana3]Manb4GlcNAcb4GlcNAc',
        'Mana2Mana6[Mana3]Mana6[Mana3]Manb4GlcNAcb4GlcNAc',
        'Mana6[Mana2Mana3]Mana6[Mana3]Manb4GlcNAcb4GlcNAc',
        'Mana6[Mana3]Mana6[Mana2Mana3]Manb4GlcNAcb4GlcNAc',
        'Mana2Mana6[Mana2Mana3]Mana6[Mana3]Manb4GlcNAcb4GlcNAc',
        'Mana2Mana6[Mana3]Mana6[Mana2Mana3]Manb4GlcNAcb4GlcNAc',
        'Mana6[Mana2Mana3]Mana6[Mana2Mana3]Manb4GlcNAcb4GlcNAc',
        'Mana2Mana6[Mana2Mana3]Mana6[Mana2Mana3]Manb4GlcNAcb4GlcNAc',
        'Mana2Mana6[Mana2Mana3]Mana6[Mana2Mana2Mana3]Manb4GlcNAcb4GlcNAc',
        'Mana6[Mana2Mana3]Mana6[Mana2Mana2Mana3]Manb4GlcNAcb4GlcNAc',
        'Mana2Mana6[Mana3]Mana6[Mana2Mana2Mana3]Manb4GlcNAcb4GlcNAc'
        ])

glycan=re.sub('Neu','Neu5Ac',glycan)
a=glycan
numbers=re.compile(r'([a|b]\d)')
a=numbers.sub(r'(\1)',a)
print(a)

