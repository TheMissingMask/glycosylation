#!/usr/bin/python

import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

def linesBetween(lines,startString,stopString):

    '''collect a list of the lines between startString and stopString'''

    outLines=[]
    copy=False

    if stopString=='':
        for l in lines:
            if l=='\n':
                next
            elif l.strip()==startString:
                copy=True
            elif copy:
                if len(l.split())<1 or l.split()[0]==';' or l.split()[0][0]==';':
                    next
                else:
                    if 1==0:#';' in l:
                        next
                    else:
                        outLines.append(l)
    else:
        for l in lines:
            if l=='\n':
               next
            elif l.strip()==startString:
                copy=True
            elif l.strip()==stopString:
              copy=False
            elif copy:
                if len(l.split())<1 or l.split()[0]==';' or l.split()[0][0]==';':
                    next
                else:
                    if 1==0:# ';' in l:
                        next
                    else:
                        outLines.append(l)

    return outLines

sns.set_style('white')

fig,ax=plt.subplots()
ax.set_aspect('equal')
G=nx.Graph()


# read the glycan itp
f=open('test.itp')
lines=f.readlines()
f.close()

# collect the residues and bonds
resiList=[]
bondList=[]
idList=[]

for l in linesBetween(lines,'[ atoms ]','[ bonds ]'):
    if l.split()[1] != 'VS1':
        next
    else:
        resiList.append(l.split()[3])
        idList.append(l.split()[0])

for l in linesBetween(lines,'[ bonds ]','[ constraints ]'):
    bondList.append((l.split()[0],l.split()[1]))

# find the branches
branchPoint=0
seen=False
stored=[]
for bond in bondList:
    firstresi=int(bond[0])
    if firstresi in stored:
        seen=firstresi
    stored.append(firstresi)
branchPoint=seen

branch0=[]
for i in range(len(idList)):
    if int(idList[i])<=branchPoint:
        branch0.append(resiList[i])

branch1Start=False
for bond in bondList:
    if int(bond[0])==branchPoint and branch1Start==False:
        b1=int(bond[1])
        branch1Start=True
    elif int(bond[0])==branchPoint and branch1Start==True:
        b2=int(bond[1])
branch1=[]
branch2=[]

b1IDs=[]
for bond in bondList[4:]:
    if str(b2) in bond:
        break
    else:
        if int(bond[0])>=b1:
            b1IDs.append(bond[0])
            b1IDs.append(bond[1])
b1IDs=list(set(b1IDs))
for i in range(len(idList)):
    if idList[i] in b1IDs:
        branch1.append(resiList[i])
b2IDs=[]
for bond in bondList[4:]:
    if bond[0] in b1IDs or bond[1] in b1IDs:
        next
    else:
        if int(bond[0])>=b2 and int(bond[1])>=b2:
            b2IDs.append(bond[0])
            b2IDs.append(bond[1])
b2IDs=list(set(b2IDs))
for i in range(len(idList)):
    if idList[i] in b2IDs:
        branch2.append(resiList[i])

#create the list of CG residues
resiTranslation={
        'GNC':'NHX',
        'MAN':'HXS',
        'GLC':'HXS',
        'FUC':'DHX',
        'GAL':'HXS'
        }

ATString=''
glycanString=''

for resi in branch0:
    glycanString+='%s-'%(resiTranslation[resi])
    ATString+='%s-'%(resi)
glycanString=glycanString[:-1]
ATString=ATString[:-1]
glycanString+='('
ATString+='('
for resi in branch1:
    glycanString+='%s-'%(resiTranslation[resi])
    ATString+='%s-'%(resi)
glycanString=glycanString[:-1]+')('
ATString=ATString[:-1]+')('
for resi in branch2:
    glycanString+='%s-'%(resiTranslation[resi])
    ATString+='%s-'%(resi)
glycanString=glycanString[:-1]+')'
ATString=ATString[:-1]+')'

#check for a fucose branch
fucBond=False
for i in idList:
    if int(i)<branchPoint:
        next
    else:
        if i in bondList[-1]:
            fucBond=bondList[-1]

if fucBond:
    glycan=glycanString[:3]+'(DHX)'+glycanString[3:]
    ATString=ATString[:3]+'(FUC)'+ATString[3:]
else:
    glycan=glycanString

print(ATString)
print(glycan)

residues=re.findall(r"[\w']+",glycan)
nodeTypes={
        'NHX':[],
        'HXS':[],
        'SIA':[],
        'DHX':[]
        }

j=0
for i in range(len(residues)):
    r=residues[i]
    if r=='DHX':
        j=1
    else:
        nodeTypes[r].append(i-j)

if fucBond:
    nodeTypes['DHX'].append(i)

if fucBond:
    branch0=glycan.split('(')[0]+'-'+glycan.split('-')[1]+'-'+glycan.split('-')[2]
    branch1=glycan.split('(')[2][:-1]
    branch2=glycan.split('(')[3][:-1]
    branch3=glycan.split('(')[1][:-1]
else:
    branch0=glycan.split('(')[0]
    branch1=glycan.split('(')[1][:-1]
    branch2=glycan.split('(')[2][:-1]

pos={}
for k in range(len(residues)):
    pos[k]=0
j=0
for i in range(len(branch0.split('-'))):
    pos[i]=(1,0+j)
    j+=1
l=0
for i in range(len(branch1.split('-'))):
    k=i+len(branch0.split('-'))
    pos[k]=(2,l+j)
    l+=1
m=0
for i in range(len(branch2.split('-'))):
    k=i+len(branch0.split('-'))+len(branch1.split('-'))
    pos[k]=(0,m+j)
    m+=1
if fucBond:
    k=0+len(branch0.split('-'))+len(branch1.split('-'))+len(branch2.split('-'))
    pos[k]=(2,0)

G.add_nodes_from(nodeTypes['NHX'],Type='NHX')
G.add_nodes_from(nodeTypes['HXS'],Type='HXS')
G.add_nodes_from(nodeTypes['SIA'],Type='SIA')
G.add_nodes_from(nodeTypes['DHX'],Type='DHX')

for n,p in pos.items():
    G.nodes[n]['pos']=p

for i in range(len(branch0.split('-'))-1):
    G.add_edge(i,i+1)
for i in range(len(branch1.split('-'))-1):
    k=i+len(branch0.split('-'))
    G.add_edge(k,k+1)
for i in range(len(branch2.split('-'))-1):
    k=i+len(branch0.split('-'))+len(branch1.split('-'))
    G.add_edge(k,k+1)
if fucBond:
    k=0+len(branch0.split('-'))+len(branch1.split('-'))+len(branch2.split('-'))
    G.add_edge(0,k)

G.add_edge(len(branch0.split('-'))-1,len(branch0.split('-')))
G.add_edge(len(branch0.split('-'))-1,len(branch0.split('-'))+len(branch1.split('-')))
    
NHX_nodes=[n for (n,ty) in nx.get_node_attributes(G,'Type').items() if ty=='NHX']
HXS_nodes=[n for (n,ty) in nx.get_node_attributes(G,'Type').items() if ty=='HXS']
SIA_nodes=[n for (n,ty) in nx.get_node_attributes(G,'Type').items() if ty=='SIA']
DHX_nodes=[n for (n,ty) in nx.get_node_attributes(G,'Type').items() if ty=='DHX']

#pos=nx.get_node_attributes(G,'pos')
nx.draw_networkx_nodes(G,pos,nodelist=NHX_nodes,node_color='xkcd:indigo',node_shape='s',label='NHX')
nx.draw_networkx_nodes(G,pos,nodelist=HXS_nodes,node_color='xkcd:aquamarine',node_shape='o',label='HXS')
nx.draw_networkx_nodes(G,pos,nodelist=SIA_nodes,node_color='xkcd:brick red',node_shape='D',label='SIA')
nx.draw_networkx_nodes(G,pos,nodelist=DHX_nodes,node_color='xkcd:grey blue',node_shape='^',label='DHX')
nx.draw_networkx_edges(G,pos,edgelist=G.edges)

minimum=np.min((ax.get_xlim(),ax.get_ylim()))
maximum=np.max((ax.get_xlim(),ax.get_ylim()))

ax.set_xlim(minimum*1.2,maximum*1.2)
ax.set_ylim(minimum*1.2,maximum*1.2)

ax.set_xticks([])
ax.set_yticks([])
ax.set_ylabel('')
ax.set_xlabel('')

plt.legend(scatterpoints=1)

plt.savefig('CG-glycan.png',dpi=600)
plt.close()

print('Gluppit the prawling strangles, there!')
