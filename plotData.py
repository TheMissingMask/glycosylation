#!/usr/bin/python

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import pandas as pd

f=open('data.dat')
lines=f.readlines()
f.close()

probabilityDict={}
for l in lines:
    cols=l.split()
    probabilityDict[cols[0]]=float(cols[1])

fig,ax=plt.subplots(figsize=(8.27,11.69),dpi=600)
sns.set_style('white')
sns.set_palette('dark')

data1=[]
for k in probabilityDict.keys():
    if probabilityDict[k]>=0.01:
        next
    elif probabilityDict[k]>=0.001:
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

# now remove all the question marked ones


editedDict={}
for k in probabilityDict.keys():
    if '?' in k:
        next
    else:
        editedDict[k]=probabilityDict[k]

fig,ax=plt.subplots(figsize=(8.27,11.69),dpi=600)
sns.set_style('white')
sns.set_palette('dark')

data1=[]
for k in editedDict.keys():
    if editedDict[k]>=0.004:
        next
    elif editedDict[k]>=0.001:
        data1.append([k,editedDict[k]])
df1=pd.DataFrame(data1,columns=['pair','probability'])
sns.barplot(data=df1,y='pair',x='probability',ax=ax,palette='GnBu_d')

ax.set_xlabel('')
ax.set_ylabel('')
ax.tick_params(axis='both', which='major', labelsize=7)
ax.tick_params(axis='both', which='minor', labelsize=7)

plt.tight_layout()
plt.savefig('edited_lower-probabilities-all.png',dpi=600)
plt.close()

fig,ax=plt.subplots(figsize=(8.27,11.69),dpi=600)
sns.set_style('white')
sns.set_palette('dark')

data1=[]
for k in editedDict.keys():
    if editedDict[k]<0.004:
        next
    else:
        data1.append([k,editedDict[k]])
df1=pd.DataFrame(data1,columns=['pair','probability'])
sns.barplot(data=df1,y='pair',x='probability',ax=ax,palette='GnBu_d')

ax.set_xlabel('')
ax.set_ylabel('')
ax.tick_params(axis='both', which='major', labelsize=7)
ax.tick_params(axis='both', which='minor', labelsize=7)

plt.tight_layout()
plt.savefig('edited_higher-probabilities-all.png',dpi=600)
plt.close()

f=open('data-edited.dat','w')
for k in editedDict.keys():
    f.write('%s   %s\n'%(k,editedDict[k]))
f.close()

print('He that would make a pun would pick a pocket')
