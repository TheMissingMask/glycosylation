#!/usr/bin/python

import re
import random
import numpy as np
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('-t','--glycanType')
parser.add_argument('-o','--outName')
args=parser.parse_args()

def extendBranch(branch):#,cappingResidues,pairProbabilities,rresProbabilities):

    tP=len(pairProbabilities.keys())/(len(pairProbabilities.keys())+len(rresProbabilities.keys()))
    i=0
    counter=0
    nres=branch.split('-')[0][:-2]

    capped=False

    while capped==False:
        for pair in pairProbabilities.keys():
            if pair.split('-')[1][1:]==nres:
                PA=pairProbabilities[pair]
                PB=rresProbabilities[nres]
                PC=PA/PB
                tP-=i
                rN=random.uniform(0,tP)
                if PC>=rN:
                    counter+=1
                    newPair=pair.split('-')[0]+'-'+pair.split('-')[1][0]
                    branch=newPair+branch
                    for cap in cappingResidues:
                        if cap in pair.split('-'):
                            capped=True
                            return branch
            addto=random.choice([0.0001,0.0005,0.001,0.0015,0.002])
            if counter>=2:
                addto=random.choice([0.005,0.01])
            i+=addto
def removePattern(lst,pattern):
    return [re.sub(pattern,'',i) for i in lst]

def tidyBonds(glycan):

    bonds=re.findall(r'\w\d-\d',glycan)

    bonds=removePattern(bonds,'1-')
    bonds=removePattern(bonds,'2-')
    return bonds
def tidyResidues(glycan):

    residues=re.split('-',str(glycan))
    residues=removePattern(residues,'a1')
    residues=removePattern(residues,'a2')
    residues=removePattern(residues,'b1')
    residues=removePattern(residues,'6\(')
    residues=removePattern(residues,'3\)')
    residues=removePattern(residues,'[0-4]')
    residues=removePattern(residues,'[6-9]')

    return residues

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


rShape={
        'Neu5Ac':'D',
        'Gal':'o',
        'Man':'o',
        'GlcNAc':'s',
        'GalNAc':'s',
        'Fuc':'^'
        }
rColour={
        'Neu5Ac':'xkcd:purple',
        'Gal':'xkcd:yellow',
        'Man':'xkcd:green',
        'GlcNAc':'xkcd:blue',
        'GalNAc':'xkcd:yellow',
        'Fuc':'xkcd:red'
        }

if args.glycanType=='oligomannose':

    possibleStructures={
            0:'Mana1-6(Mana1-3)Mana1-6(Mana1-3)Manb1-4GlcNAcb1-4GlcNAc',
            1:'Mana1-6(Mana1-3)Manb1-4GlcNAcb1-4GlcNAc',
            2:'Mana1-6(Mana1-3)Mana1-6(Mana1-2Mana1-3)Manb1-4GlcNAcb1-4GlcNAc',
            3:'Mana1-6(Mana1-2Mana1-3)Mana1-6(Mana1-3)Manb1-4GlcNAcb1-4GlcNAc',
            4:'Mana1-2Mana1-6(Mana1-3)Mana1-6(Mana1-2Mana1-3)Manb1-4GlcNAcb1-4GlcNAc',
            5:'Mana1-6(Mana1-2Mana1-3)Mana1-6(Mana1-2Mana1-3)Manb1-4GlcNAcb1-4GlcNAc',
            6:'Mana1-6(Mana1-3)Mana1-6(Mana1-2Mana1-2Mana1-3)Manb1-4GlcNAcb1-4GlcNAc',
            7:'Mana1-2Mana1-6(Mana1-2Mana1-3)Mana1-6(Mana1-2Mana1-3)Manb1-4GlcNAcb1-4GlcNAc',
            8:'Mana1-2Mana1-6(Mana1-3)Mana1-6(Mana1-2Mana1-2Mana1-3)Manb1-4GlcNAcb1-4GlcNAc',
            9:'Mana1-6(Mana1-2Mana1-3)Mana1-6(Mana1-2Mana1-2Mana1-3)Manb1-4GlcNAcb1-4GlcNAc',
            10:'Mana1-2Mana1-6(Mana1-2Mana1-3)Mana1-6(Mana1-2Mana1-2Mana1-3)Manb1-4GlcNAcb1-4GlcNAc'
            }
    n=10
    pos={
        1:{
            0:(1,1), #GlcNAc
            1:(1,2), #GlcNAc
            2:(1,3), #Man
            3:(2,4), #Man--6
            4:(2.5,5), #Man--6
            5:(1.5,5), #Man--3
            6:(0,4) #Man---3
            },
        2:{
            0:(1,1), #GlcNAc
            1:(1,2), #GlcNAc
            2:(1,3), #Man
            3:(2,4), #Man--6
            4:(2.5,5), #Man--6
            5:(1.5,5), #Man--3
            6:(0,4), #Man---3
            7:(0,5) #Man--2
            },
        3:{
            0:(1,1), #GlcNAc
            1:(1,2), #GlcNAc
            2:(1,3), #Man
            3:(2,4), #Man--6
            4:(2.5,5), #Man--6
            5:(1.5,5), #Man--3
            6:(1.5,6), #Man--2
            7:(0,4) #Man---3
            },
        4:{
            0:(1,1), #GlcNAc
            1:(1,2), #GlcNAc
            2:(1,3), #Man
            3:(2,4), #Man--6
            4:(2.5,5), #Man--6
            5:(2.5,6), #Man--2
            6:(1.5,5), #Man--3
            7:(0,4), #Man---3
            8:(0,5) #Man--2
            },
        5:{
            0:(1,1), #GlcNAc
            1:(1,2), #GlcNAc
            2:(1,3), #Man
            3:(2,4), #Man--6
            4:(2.5,5), #Man--6
            5:(1.5,5), #Man--3
            6:(1.5,6), #Man--2
            7:(0,4), #Man---3
            8:(0,5) #Man--2
            },
        6:{
            0:(1,1), #GlcNAc
            1:(1,2), #GlcNAc
            2:(1,3), #Man
            3:(2,4), #Man--6
            4:(2.5,5), #Man--6
            5:(1.5,5), #Man--3
            6:(0,4), #Man---3
            7:(0,5), #Man--2
            8:(0,6) #Man--2
            },
        7:{
            0:(1,1), #GlcNAc
            1:(1,2), #GlcNAc
            2:(1,3), #Man
            3:(2,4), #Man--6
            4:(2.5,5), #Man--6
            5:(2.5,6), #Man--2
            6:(1.5,5), #Man--3
            7:(1.5,6), #Man--2
            8:(0,4), #Man---3
            9:(0,5) #Man--2
            },
        8:{
            0:(1,1), #GlcNAc
            1:(1,2), #GlcNAc
            2:(1,3), #Man
            3:(2,4), #Man--6
            4:(2.5,5), #Man--6
            5:(2.5,6), #Man--2
            6:(1.5,5), #Man--3
            7:(0,4), #Man---3
            8:(0,5), #Man--2
            9:(0,6) #Man--2
            },
        9:{
            0:(1,1), #GlcNAc
            1:(1,2), #GlcNAc
            2:(1,3), #Man
            3:(2,4), #Man--6
            4:(2.5,5), #Man--6
            5:(1.5,5), #Man--3
            6:(1.5,6), #Man--2
            7:(0,4), #Man---3
            8:(0,5), #Man--2
            9:(0,6) #Man--2
            },
        10:{
            0:(1,1), #GlcNAc
            1:(1,2), #GlcNAc
            2:(1,3), #Man
            3:(2,4), #Man--6
            4:(2.5,5), #Man--6
            5:(2.5,6), #Man--2
            6:(1.5,5), #Man--3
            7:(1.5,6), #Man--2
            8:(0,4), #Man---3
            9:(0,5), #Man--2
            10:(0,6) #Man--2
            }
        }
  
    bondDict={
            0:['b4','b4','a6','a3','a3'],
            10:['b4','b4','a6','a2','a3','a2','a3','a2','a2']
            }
    bondPos={
        1:{
            0:(1,1.5),
            1:(1,2.5),
            2:(1.5,3.5),
            3:(2,4.5),
            4:(1.5,4.5),
            5:(0.5,3.5)
            },
        2:{
            },
        3:{
            },
        4:{
            },
        5:{
            },
        6:{
            },
        7:{
            },
        8:{
            },
        9:{
            },
        10:{
            0:(1,1.5),
            1:(1,2.5),
            2:(1,3.5),
            3:(1.25,4.5),
            4:(2,5.5),
            5:(1.5,4.5),
            6:(1.5,4.5),
            7:(0.5,3.5),
            8:(0,4.5)
            }
        }


   ###################
    resDict={
            0:['GlcNAc','GlcNAc','Man','Man','Man','Man','Man'],
            1:['GlcNAc','GlcNAc','Man','Man','Man','Man','Man','Man'],
            2:['GlcNAc','GlcNAc','Man','Man','Man','Man','Man','Man'],
            3:['GlcNAc','GlcNAc','Man','Man','Man','Man','Man','Man'],
            4:['GlcNAc','GlcNAc','Man','Man','Man','Man','Man','Man','Man'],
            5:['GlcNAc','GlcNAc','Man','Man','Man','Man','Man','Man','Man'],
            6:['GlcNAc','GlcNAc','Man','Man','Man','Man','Man','Man','Man'],
            7:['GlcNAc','GlcNAc','Man','Man','Man','Man','Man','Man','Man','Man'],
            8:['GlcNAc','GlcNAc','Man','Man','Man','Man','Man','Man','Man','Man'],
            9:['GlcNAc','GlcNAc','Man','Man','Man','Man','Man','Man','Man','Man'],
            10:['GlcNAc','GlcNAc','Man','Man','Man','Man','Man','Man','Man','Man','Man']
            }
    sns.set_style('white')

    fig,ax=plt.subplots()
    ax.set_aspect('equal')

    for r in ['GlcNAc','Man']:
        x=[]
        y=[]
        for i in range(len(resDict[n])):
            if resDict[n][i]==r:
                x.append(pos[n][i][0])
                y.append(pos[n][i][1])
        plt.scatter(x,y,marker=rShape[r],color=rColour[r])

    minimum=np.min((ax.get_xlim(),ax.get_ylim()))
    maximum=np.max((ax.get_xlim(),ax.get_ylim()))

    ax.set_xlim(minimum*1.2-1,maximum*1.2)
    ax.set_ylim(minimum*1.2,maximum*1.2)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel('')
    ax.set_xlabel('')

    fig.patch.set_visible(False)
    ax.axis('off')

    for txt in ['b4','a3','a6','a2']:
        for i in range(len(bondDict[n])):
            if bondDict[n][i]==txt:
                x=bondPos[n][i][0]
                y=bondPos[n][i][1]
                ax.annotate(txt,(x,y))
    fig.tight_layout()
#    print(newGlycan)
#    plt.savefig('%s.png'%(args.outName))
#    plt.close()
    plt.show()
        

elif args.glycanType=='hybrid':
    print('hybrid')

elif args.glycanType=='complex':
    print('complex')

    f=open('CFG-probabilities-N.dat')
    lines=f.readlines()
    f.close()

    pairProbabilities={}
    for l in lines:
        cols=l.split()
        pairProbabilities[cols[0]]=float(cols[1])


    rresProbabilities={}
    for l in lines:
        cols=l.split()
        rres=cols[0].split('-')[1][1:]
        if rres not in rresProbabilities.keys():
            rresProbabilities[rres]=0
        rresProbabilities[rres]+=float(cols[1])

    cappingResidues=[
        'Fuca1',
        'Gala1',
        'GlcNAca1',
        'Neu5Aca2'
        ]

    test='GlcNAcb1-2Mana1-6(GlcNAcb1-2Mana1-3)Manb1-4GlcNAcb1-4GlcNAc'
    branch0=test.split('(')[0]
    branch1=test.split('(')[1].split(')')[0]
    core=test.split(')')[1]


    newGlycan='%s(%s)%s'%(extendBranch(branch0),extendBranch(branch1),core)

    branch0=newGlycan.split('(')[0]
    branch1=newGlycan.split('(')[1].split(')')[0]
    core=newGlycan.split(')')[1]

# start with the N-linked core as an nx graph


    nodeNumbers=[]
    nodeTypes={
        'Neu5Ac':[],
        'Gal':[],
        'GlcNAc':[],
        'GalNAc':[],
        'Fuc':[],
        'Man':[]
        }

    i=0

    coreResidues=tidyResidues(core)[::-1]
    branch0Residues=tidyResidues(branch0[:-2])[::-1]
    branch1Residues=tidyResidues(branch1[:-2])[::-1]

    j=0
    for residues in [coreResidues,branch0Residues,branch1Residues]:
        for i in range(len(residues)):
            nodeTypes[residues[i]].append(j)
            nodeNumbers.append(j)
            j+=1

# can similarly make a list of bonds
#bondTypes={
#       'a1':[],
#       'b1':[],
#       'a2':[],
#       'b2':[],
#       'a3':[],
#       'b3':[],
#       'a4':[],
#       'b4':[],
#       'a6':[],
#       'b6':[],
#       'a8':[],
#       'b8':[]
#        }


    coreBonds=tidyBonds(core)[::-1]
    branch0Bonds=tidyBonds(branch0)[::-1]
    branch1Bonds=tidyBonds(branch1)[::-1]
#j=0
#for bonds in [coreBonds,branch0Bonds,branch1Bonds]:
#    for i in range(len(bonds)):
#        bondTypes[bonds[i]].append((j,j+1))
#        j+=1


# now we have dictionaries including identifiers for all residues and bonds from the newGlycan string -- next, we need to construct the graph


#nodeNumbers=sorted(nodeNumbers)

    branch0=newGlycan.split('(')[0]
    branch1=newGlycan.split('(')[1].split(')')[0]
    core=newGlycan.split(')')[1]

    x=1
    y=0
    pos={}
    bondPos={}

    residues=coreResidues+branch0Residues+branch1Residues
    for i in range(len(residues)):
        if i==len(core.split('-')):
            x=2
        if i==len(core.split('-'))+len(branch0.split('-'))-1:
            x=0
            y=len(core.split('-'))
        y+=1
        pos[i]=(x,y)

    x=1
    y=0.5
    for i in range(len(residues)):
        if i==len(core.split('-'))-1:
            x=1.5
        elif i>=len(core.split('-')) and i<len(core.split('-'))+len(branch0.split('-'))-2:
            x=2
        elif i==len(core.split('-'))+len(branch0.split('-'))-2:
            x=0.5
            y=len(core.split('-'))-0.5
        elif i>=len(core.split('-'))+len(branch0.split('-'))-1:
            x=0
        y+=1
        bondPos[i]=(x,y-0.2)

        bondNames=coreBonds+branch0Bonds+branch1Bonds



    sns.set_style('white')

    fig,ax=plt.subplots()
    ax.set_aspect('equal')


    for r in list(set(residues)):
        x=[]
        y=[]
        for i in nodeTypes[r]:
            x.append(pos[i][0])
            y.append(pos[i][1])
        plt.scatter(x,y,marker=rShape[r],color=rColour[r])

    minimum=np.min((ax.get_xlim(),ax.get_ylim()))
    maximum=np.max((ax.get_xlim(),ax.get_ylim()))

    ax.set_xlim(minimum*1.2-1,maximum*1.2)
    ax.set_ylim(minimum*1.2,maximum*1.2)



    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel('')
    ax.set_xlabel('')

    fig.patch.set_visible(False)
    ax.axis('off')

    for i,txt in enumerate(bondNames):
        ax.annotate(txt,bondPos[i])

    fig.tight_layout()
    print(newGlycan)
    plt.savefig('%s.png'%(args.outName))
    plt.close()
