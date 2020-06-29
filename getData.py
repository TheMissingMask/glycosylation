#!/usr/bin/python

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import pandas as pd

def tidyPair(pair):
    pair=re.sub(r'\)','',pair)
    pair=re.sub(r'\(','',pair)
    pair=re.sub(r'\]','',pair)
    pair=re.sub(r'\[','',pair)
    pair=re.sub(r'-','',pair)
    pair=re.sub(r"'",'',pair)
    pair=re.sub(r'\\n','',pair)
    pair=re.sub(r'\?','',pair)
#    if len(pair)==0:
#        return pair
#    elif re.split('\d',pair)[1]!='Rha' and pair[-1] in ['a','b']:# and re.split('\d',pair)[1]!='Rha' and len(pair)>1:
#        pair=pair[:-1]
#    elif re.split('\d',pair)[1]=='Rhaa':
#        pair=pair[:-1]
    return pair

seqList=[]
f=open('glytoucan-data')
lines=f.readlines()
f.close()
for l in lines:
    if 'WURCS' not in l or '[info]' in l:
        next
    else:
        seq=l.split(':')[1].strip()[1:-1]
        seqList.append(seq)

glycanList=[]

for seq in seqList:
    cmd='java -jar target/glycanformatconverter.jar -i WURCS -e IUPAC-Short -seq %s'%(seq)
    process=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    iupac,err=process.communicate()
    glycanList.append(str(iupac))

pairDict={}
branches=[]

totalCount=0

for glycan in glycanList:

    branches=re.findall('\(.*?\)',glycan)
    branches.append(re.sub('\(.*?\)','',glycan))

    for branch in branches:
        residues=re.split(r'[^5\d\?]',branch)
        bonds=re.findall(r'[^5\d\?]',branch)
        if residues[-1]==')':
            residues=residues[:-1]
        elif residues[0]=='(':
            residues=residues[1:]
        elif residues[0]=='[':
            residues=residues[1:]
        elif residues[0]==']':
            residues=residues[1:]
        if len(residues)==1:
            next
#        for i in range(len(residues)):
#            if residues[i]=='Neu':
#                residues[i]='Neu5Ac'
        else:
            for i in range(len(residues)-1):
                pair='%s%s%s'%(residues[i],bonds[i],residues[i+1][:-1])
                if len(pair.split(')'))>1 or len(pair.split('('))>1:
                    next
                elif 'groDman' in pair:
                    next
                elif pair[-1] in ['-','a','b',"'",'\(','\)']:
                    pair=pair[:-1]
                elif pair[0] in ["'",'\(','\)','a','b','-']:
                    pair=pair[1:]
                pair=tidyPair(pair)
                if len(pair)==0:
                    next
                elif pair not in pairDict.keys():
                    pairDict[pair]=1
                else:
                    pairDict[pair]+=1
                    totalCount+=1


# to combine and to remove

for k in pairDict.keys():
    print('%s --- %s\n'%(k,pairDict[k]))


probabilityDict={}
for k in pairDict.keys():
    probabilityDict[k]=int(pairDict[k])/totalCount

fig,ax=plt.subplots(figsize=(8.27,11.69),dpi=600)
sns.set_style('white')
sns.set_palette('dark')

data1=[]
for k in probabilityDict.keys():
    if probabilityDict[k]>=0.01:
        next
    else:
        data1.append([k,probabilityDict[k]])
df1=pd.DataFrame(data1,columns=['pair','probability'])
sns.barplot(data=df1,y='pair',x='probability',ax=ax,palette='GnBu_d')

ax.set_xlabel('')
ax.set_ylabel('')
ax.tick_params(axis='both', which='major', labelsize=7)
ax.tick_params(axis='both', which='minor', labelsize=7)

plt.tight_layout()
plt.savefig('lower-probabilities-all.png',dpi=600)
plt.close()

fig,ax=plt.subplots(figsize=(8.27,11.69),dpi=600)
sns.set_style('white')
sns.set_palette('dark')

data1=[]
for k in probabilityDict.keys():
    if probabilityDict[k]<0.01:
        next
    else:
        data1.append([k,probabilityDict[k]])
df1=pd.DataFrame(data1,columns=['pair','probability'])
sns.barplot(data=df1,y='pair',x='probability',ax=ax,palette='GnBu_d')

ax.set_xlabel('')
ax.set_ylabel('')
ax.tick_params(axis='both', which='major', labelsize=7)
ax.tick_params(axis='both', which='minor', labelsize=7)

plt.tight_layout()
plt.savefig('higher-probabilities-all.png',dpi=600)
plt.close()

f=open('data.dat','w')
for k in probabilityDict.keys():
    f.write('%s   %s\n'%(k,probabilityDict[k]))
f.close()

print('He that would make a pun would pick a pocket')
