#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import pandas as pd

def removePreNumbers(string):
    if string[0].isdigit():
        string=string[1:]
    if string[-1].isdigit():
        string=string[:-1]
    return string

def removeWhiteSpace(string):
    string=re.sub(r'\s+','',string)
    return string

def removeTermina(string):
    if string[-1] in ['a','b']:
        string=string[:-1]
    return string

def removeBrackets(string):
    if string[0] in ['(',')']:
        string=string[1:]
    if string[-1] in ['(',')']:
        string=string[:-1]
    string=re.sub(r'\(','',string)
    string=re.sub(r'\)','',string)
    return string

def getBranchPairs(branch,pairDict):

    splitBranch=branch.split('-')

    for i in range(len(splitBranch)-1):
        if len(splitBranch[i+1])<=1 or 'Sp' in splitBranch[i+1] or 'SP' in splitBranch[i+1]:
            next
        else:
            pair='%s-%s'%(splitBranch[i],splitBranch[i+1])
            pair=removePreNumbers(pair)
            pair=removeWhiteSpace(pair)
            pair=removeTermina(pair)
            pair=removeBrackets(pair)
            if pair not in pairDict.keys() and pair not in ['a1-3Man','Man-']:
                pairDict[pair]=0
            if pair not in ['a1-3Man','Man-']:
                pairDict[pair]+=1

    return pairDict


f=open('CFGArray.csv')
lines=f.readlines()
f.close()
    
glycanList=[]
for l in lines:
    glycanList.append(l.split(',')[1])
    
print('non-reducing is closer to the end e.g. Gal in LacNAc')
print('reducing is closer to the start e.g. GlcNAc in LacNAc')

pairDict={}
for glycan in glycanList:
    branchTest=glycan.split('(')
    if len(branchTest)==1:
        branch0=glycan.split('-S')[0]
        if len(branch0.split('-'))==1:
            next
        else:
            pairDict=getBranchPairs(branch0,pairDict)
#            for i in range(len(glycan.split('-'))-1):
#                pair=('%s-%s'%(branch0.split('-')[i],branch0.split('-')[i+1]))
#                if pair not in pairDict.keys():
#                    pairDict[pair]=0
    elif len(branchTest)==2:
        branch0=glycan.split('(')[0]
        branch1=glycan.split('(')[1].split(')')[0]
        branch2=glycan.split(')')[1]
        pairDict=getBranchPairs(branch0,pairDict)
        pairDict=getBranchPairs(branch1,pairDict)
        pairDict=getBranchPairs(branch2,pairDict)
        branch01=branch0.split('-')[-2]+'-'+branch0.split('-')[-1][0]+branch2.split('-')[0]
        pairDict=getBranchPairs(branch01,pairDict)
        branch02=branch1.split('-')[-2]+'-'+branch1.split('-')[-1][0]+branch2.split('-')[0]
        pairDict=getBranchPairs(branch02,pairDict)
    elif len(branchTest)==3:
        branch0=glycan.split('(')[0]
        branch1=glycan.split('(')[1].split(')')[0]
        branch2=glycan.split(')')[1].split('(')[0]
        branch3=glycan.split('(')[2].split(')')[0]
        branch4=glycan.split(')')[2]
        pairDict=getBranchPairs(branch0,pairDict)
        pairDict=getBranchPairs(branch1,pairDict)
        pairDict=getBranchPairs(branch2,pairDict)
        pairDict=getBranchPairs(branch3,pairDict)
        pairDict=getBranchPairs(branch4,pairDict)
        if len(branch0.split('-')[-1])>1:
            branch01=branch0.split('-')[-2]+'-'+branch0.split('-')[-1][0]+branch2.split('-')[0]
        else:
            branch01=branch0.split('-')[0]+'-'+branch0.split('-')[-1][0]+branch2.split('-')[0]
        pairDict=getBranchPairs(branch01,pairDict)
        if len(branch1.split('-')[-1])>1:
            branch02=branch1.split('-')[-2]+'-'+branch1.split('-')[-1][0]+branch2.split('-')[0]
        else:
            branch02=branch1.split('-')[0]+'-'+branch1.split('-')[-1][0]+branch2.split('-')[0]
        pairDict=getBranchPairs(branch02,pairDict)
        if len(branch2.split('-')[-1])>1:
            branch11=branch2.split('-')[-2]+'-'+branch2.split('-')[-1][0]+branch4.split('-')[0]
            pairDict=getBranchPairs(branch11,pairDict)
        #else:
        #    branch11=branch2.split('-')[0]+'-'+branch2.split('-')[-1][0]+branch4.split('-')[0]
        #pairDict=getBranchPairs(branch11,pairDict)
        if len(branch3.split('-')[-1])>1:
            branch12=branch3.split('-')[-2]+'-'+branch3.split('-')[-1][0]+branch4.split('-')[0]
            pairDict=getBranchPairs(branch12,pairDict)
        #else:
        #    branch12=branch3.split('-')[0]+'-'+branch3.split('-')[-1][0]+branch4.split('-')[0]
#        pairDict=getBranchPairs(branch12,pairDict)
    elif len(branchTest)==4:
        branch0=glycan.split('(')[0]
        branch1=glycan.split('(')[1].split(')')[0]
        branch2=glycan.split(')')[1].split('(')[0]
        branch3=glycan.split('(')[2].split(')')[0]
        branch4=glycan.split(')')[2]
        branch5=glycan.split('(')[3].split(')')[0]
        branch6=glycan.split(')')[3]
        pairDict=getBranchPairs(branch0,pairDict)
        pairDict=getBranchPairs(branch1,pairDict)
        pairDict=getBranchPairs(branch2,pairDict)
        pairDict=getBranchPairs(branch3,pairDict)
        pairDict=getBranchPairs(branch4,pairDict)
        pairDict=getBranchPairs(branch5,pairDict)
        pairDict=getBranchPairs(branch6,pairDict)
    elif len(branchTest)==5:
        branch0=glycan.split('(')[0]
        branch1=glycan.split('(')[1].split(')')[0]
        branch2=glycan.split(')')[1].split('(')[0]
        branch3=glycan.split('(')[2].split(')')[0]
        branch4=glycan.split(')')[2]
        branch5=glycan.split('(')[3].split(')')[0]
        branch6=glycan.split(')')[3]
        branch7=glycan.split('(')[4].split(')')[0]
        branch8=glycan.split(')')[4]
        pairDict=getBranchPairs(branch0,pairDict)
        pairDict=getBranchPairs(branch1,pairDict)
        pairDict=getBranchPairs(branch2,pairDict)
        pairDict=getBranchPairs(branch3,pairDict)
        pairDict=getBranchPairs(branch4,pairDict)
        pairDict=getBranchPairs(branch5,pairDict)
        pairDict=getBranchPairs(branch6,pairDict)
        pairDict=getBranchPairs(branch7,pairDict)
        pairDict=getBranchPairs(branch8,pairDict)
    elif len(branchTest)==6:
        branch0=glycan.split('(')[0]
        branch1=glycan.split('(')[1].split(')')[0]
        branch2=glycan.split(')')[1].split('(')[0]
        branch3=glycan.split('(')[2].split(')')[0]
        branch4=glycan.split(')')[2]
        branch5=glycan.split('(')[3].split(')')[0]
        branch6=glycan.split(')')[3]
        branch7=glycan.split('(')[4].split(')')[0]
        branch8=glycan.split(')')[4]
        branch9=glycan.split('(')[5].split(')')[0]
        branch10=glycan.split(')')[5]
        pairDict=getBranchPairs(branch0,pairDict)
        pairDict=getBranchPairs(branch1,pairDict)
        pairDict=getBranchPairs(branch2,pairDict)
        pairDict=getBranchPairs(branch3,pairDict)
        pairDict=getBranchPairs(branch4,pairDict)
        pairDict=getBranchPairs(branch5,pairDict)
        pairDict=getBranchPairs(branch6,pairDict)
        pairDict=getBranchPairs(branch7,pairDict)
        pairDict=getBranchPairs(branch8,pairDict)
        pairDict=getBranchPairs(branch9,pairDict)
        pairDict=getBranchPairs(branch10,pairDict)


toRemove=[
        'MurNAcb1-4GlcNAc',
        'KDNa2-6Gal',
        'KDNa2-3Gal',
        'Fucb1-3GlcNAc',
        'GlcNAcb1-3Fuc',
        'Gala1-3Fuc',
        'Galb1-4Fuc',
        'GlcNAcb1-4Fuc',
        'Fuca1-2Gal',
        'Fuca1-4GlcNAc',
        'Fuca1-2Man',
        'Glcb1-4Glc',
        'Glcb1-6Glc',
        'GlcAb1-6Gal',
        'GlcAb1-3GlcNAc',
        'Mana1-6Man', # already have the cores ready
        'Mana1-3Man', # as above
        'Manb1-4GlcNAc',
        'Mana1-2Man',
        'Mana1-4GlcNAc',
        'Neu5Aca2-2Man',
        'Neu5Acb2-6GalNAc',
        'Neu5Acb2-6Gal',
        'Neu5Aca2-6Man',
        'Neu5Aca2-8Neu5Ac', # need to change to allow polysialylation
        'Galb1-4Glc',
        'Glca1-4Glc',
        'Glca1-6Glc',
        'GlcNAcb1-2Man', ### here and below is particularly questionable
        'Fuca1-6Man',
        'Galb1-2Man'
        ]

toCombine=[
['SGalb1-4GlcNAc','Galb1-4GlcNAc'],
['Neu5Gca2-6Gal','Neu5Aca2-6Gal'],
['Neu5Gcb2-6Gal','Neu5Acb2-6Gal'],
['Neu5Gca2-6GalNAc','Neu5Aca2-6GalNAc'],
['Neu5Gca2-8Neu5Gc','Neu5Aca2-8Neu5Ac'],
['Neu5Gca2-8Neu5Ac','Neu5Aca2-8Neu5Ac'],
['Neu5Aca2-8Neu5Gc','Neu5Aca2-8Neu5Ac'],
['GlcNAcb1-4GlcNA','GlcNAcb1-4GlcNAc'],
['Galb1-4GlcNA','Galb1-4GlcNAc'],
['GlcNAb1-3Gal','GlcNAcb1-3Gal']
]

for c in toCombine:
    pairDict[c[1]]+=pairDict[c[0]]
    del pairDict[c[0]]

for r in toRemove:
    del pairDict[r]

for k in pairDict.keys():
    print('%s --- %s\n'%(k,pairDict[k]))


numberList=list(set(pairDict.values()))
numberDict={}
for n in numberList:
    numberDict[n]=[]
for k in pairDict.keys():
    numberDict[pairDict[k]].append(k)


# heatmap with text?
# we need to make the axes based on the residue pairs, perhaps
nrresList=[]
rrresList=[]
for k in pairDict.keys():
    nrresList.append(k.split('-')[0])
    rrresList.append(k.split('-')[1])
nrresList=list(set(nrresList))
rrresList=list(set(rrresList))
x=nrresList
y=rrresList

valueList=[]
for v in pairDict.values():
    valueList.append(int(v))
valueList=sorted(valueList)
tmpList=[]
for v in valueList:
    if v>5:
        tmpList.append(v)

dict1={}
dict2={}
dict3={}
dict4={}
for k in pairDict.keys():
    if int(pairDict[k])>=114:
        dict1[k]=int(pairDict[k])
    elif int(pairDict[k])<114 and int(pairDict[k])>=41:
        dict2[k]=int(pairDict[k])
    elif int(pairDict[k])<41 and int(pairDict[k])>=15:
        dict3[k]=int(pairDict[k])
    elif int(pairDict[k])<15 and int(pairDict[k])>5:
        dict4[k]=int(pairDict[k])


fig,ax=plt.subplots(nrows=4,ncols=1,figsize=(8,11))

sns.set_style('white')
sns.set_palette('dark')

data1=[]
for k in dict1.keys():
    data1.append([k,dict1[k]])
df1=pd.DataFrame(data1,columns=['pair','count'])
sns.barplot(data=df1,y='pair',x='count',ax=ax[0],palette='GnBu_d')
data2=[]
for k in dict2.keys():
    data2.append([k,dict2[k]])
df2=pd.DataFrame(data2,columns=['pair','count'])
sns.barplot(data=df2,y='pair',x='count',ax=ax[1],palette='GnBu_d')
data3=[]
for k in dict3.keys():
    data3.append([k,dict3[k]])
df3=pd.DataFrame(data3,columns=['pair','count'])
sns.barplot(data=df3,y='pair',x='count',ax=ax[2],palette='GnBu_d')
data4=[]
for k in dict4.keys():
    data4.append([k,dict4[k]])
df4=pd.DataFrame(data4,columns=['pair','count'])
sns.barplot(data=df4,y='pair',x='count',ax=ax[3],palette='GnBu_d')

ax[0].set_xlabel('')
ax[1].set_xlabel('')
ax[2].set_xlabel('')
ax[3].set_xlabel('')
ax[0].set_ylabel('')
ax[1].set_ylabel('')
ax[2].set_ylabel('')
ax[3].set_ylabel('')
ax[0].tick_params(axis='both', which='major', labelsize=7)
ax[0].tick_params(axis='both', which='minor', labelsize=7)
ax[1].tick_params(axis='both', which='major', labelsize=7)
ax[1].tick_params(axis='both', which='minor', labelsize=7)
ax[2].tick_params(axis='both', which='major', labelsize=7)
ax[2].tick_params(axis='both', which='minor', labelsize=7)
ax[3].tick_params(axis='both', which='major', labelsize=7)
ax[3].tick_params(axis='both', which='minor', labelsize=7)

plt.tight_layout()

plt.savefig('CFG-data_over5.png',dpi=600)
plt.close()

xList=range(1,6)
yList=[]
for k in xList:
    yList.append(len(numberDict[k]))
data=[]
for i in range(len(xList)):
    data.append([int(xList[i]),yList[i]])
df=pd.DataFrame(data,columns=['count','pairs'])
# plot bar chart of the ones with up to five occurences and label the bars
fig,ax=plt.subplots(figsize=(6,8))
sns.set_style('white')
sns.barplot(data=df,x='count',y='pairs',palette='GnBu_d')
ax.set_xlabel('number of occurrences')
ax.set_ylabel('number of pairs')
plt.savefig('CFG-data5less.png',dpi=600)
plt.close()

xList=[]
yList=[]
for k in sorted(numberDict.keys()):
    xList.append(k)
    yList.append(numberDict[k])
f=open('numberData.txt','w')
for i in range(len(xList)):
    string=''
    for unit in yList[i]:
        string+='%s '%(unit)
    f.write('%s   %s\n'%(xList[i],string))
f.close()

f=open('numberData.txt')
lines=f.readlines()
f.close()
pairs=[]
for l in lines:
    cols=l.split()
    for c in cols[1:]:
        pairs.append(c)
nPairs=0
for l in lines:
    cols=l.split()
    count=int(cols[0])*len(cols[1:])
    nPairs+=count

probabilityDict={}
for l in lines:
    cols=l.split()
    p=int(cols[0])/nPairs
    for pair in cols[1:]:
        probabilityDict[pair]=p

f=open('CFG-probabilities-N.dat','w')
for k in probabilityDict.keys():
    f.write('%s   %s\n'%(k,probabilityDict[k]))
f.close()

print('He that would make a pun would pick a pocket')
