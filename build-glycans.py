#!/usr/bin/python

import argparse
import pandas as pd
import numpy as np
import re
import random
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import random

# functions #

def linesFromFile(fName):

    '''just get the lines from a file and return them as a list'''

    f=open(fName,'r')
    lines=f.readlines()
    f.close()

    return lines

def linesBetween(lines,startString,stopString):

    '''collect a list of the lines between startString and stopString'''

    outLines=[]
    copy=False

    if stopString=='':
        iterLines=iter(lines)
        for l in iterLines:
            if l=='\n':
                next(iterLines)
            elif l.strip()==startString:
                copy=True
            elif copy:
                if len(l.split())<1 or l.split()[0]==';':
                    next(iterLines)
                else:
                    outLines.append(l)
    else:
        iterLines=iter(lines)
        for l in iterLines:
            if l=='\n':
                next(iterLines)
            elif l.strip()==startString:
                copy=True
            elif l.strip()==stopString:
                copy=False
            elif copy:
                if len(l.split())<1 or l.split()[0]==';':
                    next(iterLines)
                else:
                    outLines.append(l)

    return outLines


def updateITP(itpList,bondDict,angleDict,mapDict,bond0,VS0,pair,newID,newResnum,itpLines):

    '''update the itp file with the latest glycan residue'''

    # first, let us work out what is going on
    rRes=pair.split('-')[0]
    nrRes=pair.split('-')[2]
    bond1=pair.split('-')[1]
    bondAngle=bond1

    bondPair='%s-%s'%(mapDict[rRes],mapDict[nrRes])

    newLines=itpLines

    # atoms
    newAtoms=linesBetween(newLines,'[ atoms ]','[ bonds ]')
    tmp=[]
    for l in newAtoms:
        cols=l.split()
        tmp.append('%s  %s  %s  %s  %s  %s  %s\n'%(newID[cols[0]],cols[1],newResnum,cols[3],cols[4],newID[cols[0]],cols[6]))
    VS1=int(newID[cols[0]])
    VS0=int(VS0)
    newAtoms=tmp

    # bonds
    newBonds=linesBetween(newLines,'[ bonds ]','[ constraints ]')
    tmp=[]
    for l in newBonds:
        cols=l.split()
        tmp.append('%s  %s  %s  %s  %s\n'%(newID[cols[0]],newID[cols[1]],cols[2],cols[3],cols[4]))
    newBonds=tmp
    if bondPair in bondDict.keys():
        newBonds.append('%s  %s  2  %s  %s\n'%(VS1,VS0,bondDict[bondPair][0],bondDict[bondPair][1]))
    else:
        newBonds.append('%s  %s  2  %s  %s\n'%(VS1,VS0,bondDict[None][0],bondDict[None][1]))
   
    # constraints
    newConstraints=linesBetween(newLines,'[ constraints ]','[ angles ]')
    tmp=[]
    for l in newConstraints:
        cols=l.split()
        tmp.append('%s  %s  1  %s\n'%(newID[cols[0]],newID[cols[1]],cols[3]))
    newConstraints=tmp

    # angles
    newAngles=linesBetween(newLines,'[ angles ]','[ dihedrals ]')
    tmp=[]
    for l in newAngles:
        cols=l.split()
        tmp.append('%s  %s  %s  2  %s  %s\n'%(newID[cols[0]],newID[cols[1]],newID[cols[2]],cols[4],cols[5]))
    newAngles=tmp
    ### here we add the angleDict information -- need the ring ids and the virtual site ids of the residues involved
    ### virtual sites are VS1 and VS0
    ### we can figure out the requirements for the new monosaccharide using the atoms i.e.
    countingAtoms0=linesBetween(
            linesFromFile('%s.itp'%(nrRes)),
            '[ atoms ]',
            '[ bonds ]')[-1].split()[0]
    ### as for the previous monosaccharide, we need to consider something else entirely
    countingAtoms1=linesBetween(
            linesFromFile('%s.itp'%(rRes)),
            '[ atoms ]',
            '[ bonds ]')[-1].split()[0]
    ### now we can make the angles
    if bondAngle in angleDict.keys():
        newAngles.append('%s  %s  %s  2  %s  %s\n'%(int(VS0)-(int(countingAtoms0))+1,VS0,VS1,angleDict[bondAngle][0][0],angleDict[bondAngle][0][1]))
        newAngles.append('%s  %s  %s  2  %s  %s\n'%(VS0-(int(countingAtoms0))+2,VS0,VS1,angleDict[bondAngle][1][0],angleDict[bondAngle][1][1]))
        newAngles.append('%s  %s  %s  2  %s  %s\n'%(VS0-(int(countingAtoms0))+3,VS0,VS1,angleDict[bondAngle][2][0],angleDict[bondAngle][2][1]))
        newAngles.append('%s  %s  %s  2  %s  %s\n'%(VS0,VS1,VS1-(int(countingAtoms1))+1,angleDict[bondAngle][3][0],angleDict[bondAngle][3][1]))
        newAngles.append('%s  %s  %s  2  %s  %s\n'%(VS0,VS1,VS1-(int(countingAtoms1))+2,angleDict[bondAngle][4][0],angleDict[bondAngle][4][1]))
        newAngles.append('%s  %s  %s  2  %s  %s\n'%(VS0,VS1,VS1-(int(countingAtoms1))+3,angleDict[bondAngle][5][0],angleDict[bondAngle][5][1]))

    # dihedrals -- at the moment we have none of these
    newDihedrals=linesBetween(newLines,'[ dihedrals ]','[ exclusions ]')

    # exclusions -- at the moment we have none of these
    newExclusions=linesBetween(newLines,'[ exclusions ]','[ virtual_sitesn ]')

    newVirtuals=linesBetween(newLines,'[ virtual_sitesn ]','')
    tmp=[]
    for l in newVirtuals:
        cols=l.split()
        tmp.append('%s  1  %s  %s  %s\n'%(VS1,newID['1'],newID['2'],newID['3']))
    newVirtuals=tmp

    itpList[0]=itertools.chain(itpList[0],newAtoms)
    itpList[1]=itertools.chain(itpList[1],newBonds)
    itpList[2]=itertools.chain(itpList[2],newConstraints)
    itpList[3]=itertools.chain(itpList[3],newAngles)
    itpList[4]=itertools.chain(itpList[4],newDihedrals)
    itpList[5]=itertools.chain(itpList[5],newExclusions)
    itpList[6]=itertools.chain(itpList[6],newVirtuals)

    return itpList,VS1,VS1

def loadPDB(filename):

    '''get the PDB coordinates and complete lines'''

    lines=linesFromFile(filename)

    coords=[]
    for l in lines:
        if l.split()[0]=='ATOM':
            coords.append(l)

    lineDict={}
    coordDict={}

    for c in coords:
        cols=c.split()
        lineDict[int(cols[4])]=[]
        coordDict[int(cols[4])]=[]

    for c in coords:
        cols=c.split()
        lineDict[int(cols[4])].append(c)
        coordDict[int(cols[4])].append([float(x) for x in cols[5:8]])

    return lineDict,coordDict

def updatePDB(glycanCoords,rRes,resnum0,maxResnum,maxID,pdbLines):

    '''update the coordinates of the glycan'''

    pos1=np.asarray(glycanCoords[1][-1])
    pos2=np.asarray(glycanCoords[2][-1])
    #if len(glycanCoords.keys())>=3:
    #    pos3=np.asarray(glycanCoords[3][-1])
    #else:
    #    pos3=np.asarray(pos2+(pos2-pos1))
    pos3=np.asarray(pos2+(pos2-pos1))

    translateList=np.asarray([
        np.subtract(pos2,pos1),
        np.subtract(pos3,pos2)])
    translate1=np.asarray([np.mean(translateList[:,0]),np.mean(translateList[:,1]),np.mean(translateList[:,2])])
    position=np.asarray(glycanCoords[int(resnum0)][-1])

    resLines=linesFromFile('%s.pdb'%(rRes))
    VS=np.asarray([float(x) for x in resLines[-1].split()[5:8]])
    if rRes=='Fuc':
        translate1=np.asarray([translate1[1],translate1[2],translate1[0]])
    position=np.add(translate1,position)
    translate2=np.subtract(position,VS)

    pdbLines[maxResnum+1]=[]
    glycanCoords[maxResnum+1]=[]

    for l in resLines:
        maxID+=1
        cols=l.split()
        coords=np.asarray([float(x) for x in cols[5:8]])
        newCoords=np.add(translate2,coords)
        pdbLines[maxResnum+1].append('ATOM     %3d  %s %s    %2s     %.3f %.3f %.3f  1.00  1.00\n'%(maxID,cols[2],cols[3],maxResnum+1,newCoords[0],newCoords[1],newCoords[2]))
        glycanCoords[maxResnum+1].append([newCoords[0],newCoords[1],newCoords[2]])

    return pdbLines,glycanCoords

def theMaestro(itpList,pair,bond0,VS0,maxResnum,maxID,resnum0,pdbLines,glycanCoords):

    '''orchestrates much of the rest of the glycosylation script'''

    print('I lost my shoe')

    rRes=pair.split('-')[0]
    bond1=pair.split('-')[1]
    nrRes=pair.split('-')[2]

    # itp stuff
    itpLines=linesFromFile('%s.itp'%(rRes)) # changed nrRes to rRes, because I think that makes sense???
    newResnum=int(maxResnum)+1
    newID={}
    for a in linesBetween(itpLines,'[ atoms ]','[ bonds ]'):
            newID[a.split()[0]]=int(a.split()[0])+maxID
    itpList,maxID,VS0=updateITP(itpList,bondDict,angleDict,mapDict,bond0,VS0,pair,newID,newResnum,itpLines)

    # pdb stuff
    pdbLines,glycanCoords=updatePDB(glycanCoords,rRes,resnum0,maxResnum,maxID,pdbLines)

    maxResnum=newResnum

    return itpList,maxResnum,maxID,pdbLines,glycanCoords,VS0


# set up the data dictionaries

bondDict={
        'GMY-GMY':(0.57,20),
        'GMY-GNC':(0.52,20),
        'GNC-GMY':(0.52,20),
        'GNC-GNC':(0.59,20),
        'GNC-GCA':(0.54,20),
        'GCA-GNC':(0.54,20),
        'GMY-GCA':(0.54,20),
        'GCA-GMY':(0.54,20),
        'GMY-NMC':(0.50,20),
        'NMC-GMY':(0.50,20),
        'GMY-FUC':(0.56,20),
        'FUC-GMY':(0.56,20),
        'NMC-NMC':(0.64,20),
        'GNC-FUC':(0.53,20),
        'FUC-GNC':(0.53,20),
        None:(0.55,20)
        }

angleDict={
        '2':[(13,10),(82,10),(147,10),(58,10),(155,10),(74,10)],
        '3':[(89,10),(25,10),(140,10),(56,10),(142,10),(78,10)],
        '4':[(139,10),(46,10),(96,10),(46,10),(140,10),(81,10)],
        '6':[(127,10),(126,10),(21,10),(47,10),(131,10),(82,10)],
        '8':[(90,10),(148,10),(44,10),(51,10),(23,10),(43,10)]
        }

mapDict={
        'GlcA':'GCA',
        'IdoA':'GCA',
        'Glc':'GMY',
        'Man':'GMY',
        'Gal':'GMY',
        'All':'GMY',
        'Alt':'GMY',
        'GlcNAc':'GNC',
        'GalNAc':'GNC',
        'GlcN':'GCN',
        'Neu5Ac':'NMC',
        'Neu5Gc':'NMC',
        'Fuc':'FUC',
        'Rha':'FUC',
        'Qui':'FUC',
        'Xyl':'XYL',
        'Gal3S':'GMY3S',
        'GalNAc3S':'GNC3S',
        'GlcA1S':'GCA1S',
        'IdoA1S':'GCA1S',
        'GlcNAc3S':'GNC3S',
        'GlcNAc2S':'GNC2S',
        'GalNAc2S':'GNC2S',
        'GalNAc4S':'GNC4S',
        'GlcNAc4S':'GNC4S'
        }

# get the probabilities

df=pd.read_csv('CFGArray.csv',header=None)
glycans=list(df[1])

glycanPairs=[]
for glycan in glycans:
    for k in re.findall(r'\(?\)?\w{4,}\(?\)?\-\)?\(?\w{4,}\(?\)?','%s'%(glycan)):
        if 'Sp' in k or 'MurNAc' in k: # bacterial sugars will need to go in a separate script
            next
        else:
            if k[0]=='(' or k[0]==')':
                glycanPairs.append(k[1:])
            else:
                glycanPairs.append(k)

pairDict={}
for glycan in glycans:
    for k in re.findall(r'\d(\w{4,})1-(\d)(\w{4,})\d?','%s'%(glycan)):
        if k[2][-1] not in ['a','b','1','2','3','4','5','6',1,2,3,4,5,6]:
            pair='%s-%s-%s'%(k[0][:-1],k[1],k[2])
        elif k[2][-1] in ['a','b']:
            pair='%s-%s-%s'%(k[0][:-1],k[1],k[2][:-1])
        else:
            pair='%s-%s-%s'%(k[0][:-1],k[1],k[2][:-2])
        pairDict[pair]=0
    for k in re.findall(r'(Neu\d\w{3,})2-(\d)(\w{4,})\d?','%s'%(glycan)):
        if k[2][-1] not in ['a','b','1','2','3','4','5','6',1,2,3,4,5,6]:
            pair='%s-%s-%s'%(k[0][:-1],k[1],k[2])
        elif k[2][-1] in ['a','b']:
            pair='%s-%s-%s'%(k[0][:-1],k[1],k[2][:-1])
        else:
            pair='%s-%s-%s'%(k[0][:-1],k[1],k[2][:-2])
        pairDict[pair]=0


for glycan in glycans:
    for k in re.findall(r'\d(\w{4,})1-(\d)(\w{4,})\d?','%s'%(glycan)):
        if k[2][-1] not in ['a','b','1','2','3','4','5','6',1,2,3,4,5,6]:
            pair='%s-%s-%s'%(k[0][:-1],k[1],k[2])
        elif k[2][-1] in ['a','b']:
            pair='%s-%s-%s'%(k[0][:-1],k[1],k[2][:-1])
        else:
            pair='%s-%s-%s'%(k[0][:-1],k[1],k[2][:-2])
        pairDict[pair]+=1
    for k in re.findall(r'(Neu\d\w{3,})2-(\d)(\w{4,})\d?','%s'%(glycan)):
        if k[2][-1] not in ['a','b','1','2','3','4','5','6',1,2,3,4,5,6]:
            pair='%s-%s-%s'%(k[0][:-1],k[1],k[2])
        elif k[2][-1] in ['a','b']:
            pair='%s-%s-%s'%(k[0][:-1],k[1],k[2][:-1])
        else:
            pair='%s-%s-%s'%(k[0][:-1],k[1],k[2][:-2])
        pairDict[pair]+=1



# get user input #

parser=argparse.ArgumentParser()
parser.add_argument('-t','--glycanType',default='N',help='Type of glycosylation (N-linked [N], mucin-type [O], glycosaminoglycan [GAG], monosaccharide [mono])')
parser.add_argument('--cpx',action='store_true',default=False)
parser.add_argument('--hmn',action='store_true',default=False)
parser.add_argument('--mucin',action='store_true',default=False)
parser.add_argument('--glycosaminoglycan',default=None)
parser.add_argument('--mono',default=None)
parser.add_argument('-o','--outName',default='glycan')
args=parser.parse_args()


# set everything up #


if args.glycanType=='N': # all the following only applies to N-linked glycosylation

    core='man3glcnac2'
    coreLines=linesFromFile('%s.itp'%(core))
    itpList=[
            linesBetween(coreLines,'[ atoms ]','[ bonds ]'),
            linesBetween(coreLines,'[ bonds ]','[ constraints ]'),
            linesBetween(coreLines,'[ constraints ]','[ angles ]'),
            linesBetween(coreLines,'[ angles ]','[ dihedrals ]'),
            linesBetween(coreLines,'[ dihedrals ]','[ exclusions ]'),
            linesBetween(coreLines,'[ exclusions ]','[ virtual_sitesn ]'),
            linesBetween(coreLines,'[ virtual_sitesn ]','')
            ]

    hmn=args.hmn
    cpx=args.cpx

    prSet=pairDict.keys()
    allProbabilities={}
    for pair in prSet:
        if hmn==False and pair.split('-')[0]=='Man': # if the residue to be added in mannose
            print('was trying to add mannose, but not going to now')
            next
        elif cpx==False and pair.split('-')[0]!='Man':
            next
        else:
            allProbabilities[pair]=int(pairDict[pair])/float(sum(list(pairDict.values())))

    nonreducing,reducing=[],[]
    for pair in prSet:
        nonreducing.append(pair.split('-')[2])
        reducing.append(pair.split('-')[0])
    nrSet=list(set(nonreducing))
    rSet=list(set(reducing))

    nrProbabilities={}
    for mono in nrSet:
        nrProbabilities[mono]=nonreducing.count(mono)/float(len(nonreducing))

    rProbabilities={}
    for mono in rSet:
        rProbabilities[mono]=reducing.count(mono)/float(len(reducing))

    # get the branches ready

    if args.glycanType=='N':
        branches=['Man-3-Man','Man-6-Man']

    maxIterations=2

    glycanLines,glycanCoords=loadPDB('%s.pdb'%(core))
    pdbLines=glycanLines

    maxID=int(glycanLines[len(glycanLines.keys())][-1].split()[1]) # this is the current highest atom ID from the core
    maxResnum=int(glycanLines[len(glycanLines.keys())][-1].split()[4]) # this is the current highest residue number from the core

    branch1=None # by default, we will have no further branching
    resnum1=4

    # do the iterations #

    for branch in branches: # loop through the main branches

        print('everything is connected')

        # get a degree of variation regarding branch lengths

        if random.random()<0.3:
            maxIterations-=1
        elif random.random()<0.2:
            maxIterations+=1

        # get a bunch of useful info about the branch
        if branch=='Man-3-Man':
            branchCoords=glycanCoords[3][-1] # coordinates of the VS of the final residue in the branch
        elif branch=='Man-6-Man':
            branchCoords=glycanCoords[4][-1]
        else:
            branchCoords=glycanCoords[1][-1]

        ID0=glycanLines[resnum1][-1].split()[1] # atom ID of the VS of the final residue in the branch
        resnum0=resnum1
        if resnum0>12: # don't want the glycans to have too many residues
            break
        branch0=[branch.split('-')[0],resnum0,ID0,branchCoords] # list containing the above information, for convenience I suppose

        # this will allow us to do hybrid glycosylation
#        if branch==branches[1] and hmn==True:
#            cpx=False

        # this is a check to make sure we get no more mannose residues within a complex branch
        if cpx==True:
            tmp=[]
            for r in rSet:
                if r=='Man':
                    next
                else:
                    tmp.append(r)
            rSet=tmp

        # get going with the glycosylation now, methinks
        counter=0 # use this to monitor how many rounds of glycosylation we have done

        nrRes=branch.split('-')[2] # the nonreducing residue i.e. that already added
        rRes=branch.split('-')[0] # the reducing residue i.e. that being added
        bond1=branch.split('-')[1] # the bond to be created -- just a placeholder for now, so to speak
        bond0=branch.split('-')[1] # the last bond added
        if branch==branches[0] and args.glycanType=='N':
            VS0=18 # for now
        elif branch!=branches[0] and args.glycanType=='N':
            VS0=22 # for now
        else:
            VS0=maxID

        while counter<maxIterations: # keep extending the branch until we reach the maximum number of allowed iterations for that branch

            for rRes in rSet: # for each possible next added residue

                for b in [0,1,2,3,4,5,6,7,8,9]: # and each possible bond to the added residue

                    pair='%s-%s-%s'%(rRes,b,nrRes) # what is the combination we end up with
                    if pair in allProbabilities.keys(): # and is that a combination we know about?

                        probability=float(allProbabilities[pair])/nrProbabilities[nrRes] # the probability of encountering that pair of residues linked by that bond

                        if random.random()<probability: # now, if we generate a random number that is less than the probability, we will add rRes into the branch

                            print('the shark kitten attacks') # this is a very helpful thing to be outputting...

                            itpList,maxResnum,maxID,pdbLines,glycanCoords,VS0=theMaestro(itpList,pair,bond0,VS0,maxResnum,maxID,resnum0,pdbLines,glycanCoords) # call the maestro

                            counter+=1 # update the counter

                            nrRes=rRes # now we update the nrRes, so that we can look for the next residue to be added
                            bond0=bond1 # the new bond becomes the previous one
                            resnum0=int(resnum0)+1 # and this
                        
        resnum1+=counter
        maxIterations=3

    if cpx==True and random.random()<0.5: # do this out of the loop, as only want it to happen once

        itpList,maxResnum,maxID,pdbLines,glycanCoords,VSn=theMaestro(itpList,'Fuc-3-%s'%(nrRes),bond0,5,maxResnum,maxID,1,pdbLines,glycanCoords)

elif args.glycanType=='O':
    core=random.choice(['core1','core2','core3','core4','core5','core6','core7','core8'])
  
    coreLines=linesFromFile('%s.itp'%(core))
    itpList=[
            linesBetween(coreLines,'[ atoms ]','[ bonds ]'),
            linesBetween(coreLines,'[ bonds ]','[ constraints ]'),
            linesBetween(coreLines,'[ constraints ]','[ angles ]'),
            linesBetween(coreLines,'[ angles ]','[ dihedrals ]'),
            linesBetween(coreLines,'[ dihedrals ]','[ exclusions ]'),
            linesBetween(coreLines,'[ exclusions ]','[ virtual_sitesn ]'),
            linesBetween(coreLines,'[ virtual_sitesn ]','')
            ]

    glycanLines,glycanCoords=loadPDB('%s.pdb'%(core))
    pdbLines=glycanLines
    maxID=int(glycanLines[len(glycanLines.keys())][-1].split()[1]) # this is the current highest atom ID from the core
    maxResnum=int(glycanLines[len(glycanLines.keys())][-1].split()[4]) # this is the current highest residue number from the core

    resnum1=maxResnum

    if args.glycanType in ['O1','O8']:
        branches=['Gal-3-GalNAc']
    elif args.glycanType in ['O2','O4']:
        branches=['Gal-3-GalNAc','GlcNAc-6-GalNAc']
    elif args.glycanType in ['O3','O5']:
        branches=['GlcNAc-3-GalNAc']
    elif args.glycanType in ['O6','O7']:
        branches=['GlcNAc-6-GalNAc']

    for branch in branches:

        counter=0

        if branch in ['Gal-3-GalNAc','GlcNAc-3-GalNAc']:
            branchCoords=glycanCoords[1][-1] # coordinates of the VS of the final residue in the branch
        elif branch in ['GlcNAc-6-GalNAc','GlcNAc-6-GlcNAc','GlcNAc-3-GalNAc']:
            branchCoords=glycanCoords[2][-1]
        else:
            branchCoords=glycanCoords[1][-1] # coordinates of the VS of the final residue in the branch

        # at least at the moment, we shall limit the potential residues that can be added for O-linked
        addingOptions={
            0:['Gal','Neu5Ac'],
            1:['GlcNAc','GalNAc','Neu5Ac'],
            2:['Neu5Ac']
            }

        ID0=glycanLines[resnum1][-1].split()[1] # atom ID of the VS of the final residue in the branch
        resnum0=resnum1

        branch0=[branch.split('-')[0],resnum0,ID0,branchCoords] # list containing the above information, for convenience I suppose

        nrRes=branch.split('-')[2] # the nonreducing residue i.e. that already added
        rRes=branch.split('-')[0] # the reducing residue i.e. that being added
        bond1=branch.split('-')[1] # the bond to be created -- just a placeholder for now, so to speak
        bond0=branch.split('-')[1] # the last bond added
        if branch==branches[0]:
            VS0=6 # for now
        elif branch!=branches[0]:
            VS0=maxID # for now

        while counter<3:
            # extend the core
            # somehow
            # maybe
            # i don't know -_-

            rRes=random.choice(addingOptions[counter])
            itpList,maxResnum,maxID,pdbLines,glycanCoords,VS0=theMaestro(itpList,pair,bond0,VS0,maxResnum,maxID,resnum0,pdbLines,glycanCoords) # call the maestro
            bond1=3 ### not always, but will do for the moment#######################################################################################################!!!!!!!!! NEED TO CHANGE
            counter+=1 # don't forget to do this or you shall be stuck in an infinite loop
            nrRes=rRes # now we update the nrRes, so that we can look for the next residue to be added
            bond0=bond1 # the new bond becomes the previous one
            resnum0=int(resnum0)+1 # and this
            if rRes=='Neu5Ac': # we do not extend beyond sialic acid
                print('just added sialic acid and therefore I exit (hopefully not pursued by a bear)')
                break
                        
        resnum1+=counter
        maxIterations=2


##### glycosaminoglycans #### the pair is the issue -- needs to be e.g. Gal-3-GalNAc
elif args.glycanType=='GAG':
   
    if args.glycosaminoglycan in ['HA','hyaluronan','hyaluronic acid','ha']:
        core='4-GlcA-3-GlcNAc'
    else:
        core='3-Gal-3-Gal-4-Xyl'

    coreLines=linesFromFile('%s.itp'%(core))
    itpList=[
            linesBetween(coreLines,'[ atoms ]','[ bonds ]'),
            linesBetween(coreLines,'[ bonds ]','[ constraints ]'),
            linesBetween(coreLines,'[ constraints ]','[ angles ]'),
            linesBetween(coreLines,'[ angles ]','[ dihedrals ]'),
            linesBetween(coreLines,'[ dihedrals ]','[ exclusions ]'),
            linesBetween(coreLines,'[ exclusions ]','[ virtual_sitesn ]'),
            linesBetween(coreLines,'[ virtual_sitesn ]','')
            ]

    glycanLines,glycanCoords=loadPDB('%s.pdb'%(core))
    pdbLines=glycanLines
    maxID=int(glycanLines[len(glycanLines.keys())][-1].split()[1]) # this is the current highest atom ID from the core
    maxResnum=int(glycanLines[len(glycanLines.keys())][-1].split()[4]) # this is the current highest residue number from the core

    resnum1=maxResnum

    counter=0

    branchCoords=glycanCoords[max(glycanCoords.keys())][-1] # coordinates of the VS of the final residue in the branch

    if args.glycosaminoglycan in ['heparan sulphate','heparin','hs']:
        gag='hs'
    elif args.glycosaminoglycan in ['keratan sulphate','ks']:
        gag='ks'
    elif args.glycosaminoglycan in ['chondroitin sulphate','cs']:
        gag='cs'
    elif args.glycosaminoglycan in ['dermatan sulphate','ds']:
        gag='ds'
    elif args.glycosaminoglycan in ['hyaluronan','hyaluronic acid','ha']:
        gag='ha'
    else:
        raise Exception('glycosaminoglycan not recognised')

    addingOptions={
            'ha':{0:['GlcA1S','GlcA'],
                1:['GlcNAc2S','GlcNAc3S','GlcNAc4S','GlcNAc'],
                'bonds':[3,4] #########################
                },
            'cs':{0:['GlcA1S','GlcA'],
                1:['GalNAc2S','GalNAc3S','GalNAc'],
                'bonds':[3,4] ###########################
                },
            'ds':{0:['IdoA1S','IdoA'],
                1:['GalNAc2S','GalNAc3S','GalNAc'],
                'bonds':[3,4] ###########################
                },
            'ks':{0:['Gal3S','Gal'],
                1:['GlcNAc3S','GlcNAc'],
                'bonds':['4','3']#########################
                },
            'hs':{0:['GlcA1S','GlcA'],
                1:['GlcNAc2S','GlcNAc3S','GlcNAc4S','GlcNAc'],
                'bonds':[4,4]
                }
        }
        ############################### need to check the bond order -- look at alternative outcomes and compare with atomistic models
    ID0=maxID
    resnum0=resnum1

    branch=core
    branch0=[branch.split('-')[1],resnum0,ID0,branchCoords] # list containing the above information, for convenience I suppose

    nrRes=branch.split('-')[3] # the nonreducing residue i.e. that already added
    rRes=branch.split('-')[1] # the reducing residue i.e. that being added
    bond1=branch.split('-')[0] # the bond to be created -- just a placeholder for now, so to speak
    bond0=branch.split('-')[0] # the last bond added
    VS0=maxID

    gagLength=random.choice(range(60,100))
    alt=0
    while counter<gagLength:
        rRes=random.choice(addingOptions[gag][alt])
        if alt==0:
            pair='%s-%s-%s'%(rRes,addingOptions[gag]['bonds'][0],nrRes)
        elif alt==1:
            pair='%s-%s-%s'%(rRes,addingOptions[gag]['bonds'][1],nrRes)
        itpList,maxResnum,maxID,pdbLines,glycanCoords,VS0=theMaestro(itpList,pair,bond0,VS0,maxResnum,maxID,resnum0,pdbLines,glycanCoords) # call the maestro
        bond1=addingOptions[gag]['bonds'][1]
        counter+=1 # don't forget to do this or you shall be stuck in an infinite loop
        nrRes=rRes # now we update the nrRes, so that we can look for the next residue to be added
        bond0=bond1 # the new bond becomes the previous one
        resnum0=int(resnum0)+1 # and this

#        f=open('%s.pdb'%(counter),'w')
#        for k in pdbLines.keys():
#            for l in pdbLines[k]:
#                f.write(l)
#        f.close()

        if alt==0:
            alt=1
        elif alt==1:
            alt=0

elif args.glycanType=='mono':

    if args.mono==False:
        print('You need to define a monosaccharide if you want to attach one')
    else:
        if args.mono in ['Glucose','glucose','Glc']:
            mono='Glc'
        elif args.mono in ['Mannose','mannose','Man']:
            mono='Man'
        elif args.mono in ['Galactose','galactose','Gal']:
            mono='Gal'
        elif args.mono in ['GlcNAc','N-acetylglucosamine']:
            mono='GlcNAc'
        elif args.mono in ['GalNAc','N-acetylgalactosamine']:
            mono='GalNAc'

    coreLines=linesFromFile('%s.itp'%(mono))
    itpList=[
            linesBetween(coreLines,'[ atoms ]','[ bonds ]'),
            linesBetween(coreLines,'[ bonds ]','[ constraints ]'),
            linesBetween(coreLines,'[ constraints ]','[ angles ]'),
            linesBetween(coreLines,'[ angles ]','[ dihedrals ]'),
            linesBetween(coreLines,'[ dihedrals ]','[ exclusions ]'),
            linesBetween(coreLines,'[ exclusions ]','[ virtual_sitesn ]'),
            linesBetween(coreLines,'[ virtual_sitesn ]','')
            ]

    glycanLines,glycanCoords=loadPDB('%s.pdb'%(mono))
    pdbLines=glycanLines
    
else:
    print('core not found -- seek help')
    
# yield the output #

# first the itp

headerList=[
        '[ atoms ]\n',
        '[ bonds ]\n',
        '[ constraints ]\n',
        '[ angles ]\n',
        '[ dihedrals ]\n',
        '[ exclusions ]\n',
        '[ virtual_sitesn ]\n'
        ]
f=open('%s.itp'%(args.outName),'w')
f.write('[ moleculetype ]\n %s   3\n\n'%(args.outName))
for i in range(len(itpList)):
    f.write(headerList[i])
    for l in itpList[i]:
        f.write(l)
    f.write('\n\n')
f.close()

# and now the pdb
f=open('%s.pdb'%(args.outName),'w')
for k in pdbLines.keys():
    for l in pdbLines[k]:
        f.write(l)
f.close()

print('The rolling stone gets the worm')
