#!/usr/bin/python

### need to fix the bilayer thing ###

# imports 

import subprocess
import argparse
import scipy
import numpy as np
import re
import os
import random
import itertools
import MDAnalysis
from MDAnalysis.analysis.leaflet import LeafletFinder
from scipy.spatial.distance import cdist


# define stuff
def unit(v):
    x=v[0]
    y=v[1]
    z=v[2]

    if x!=0:
        x=float(x)/np.linalg.norm(x)
    if y!=0:
        y=float(y)/np.linalg.norm(y)
    if z!=0:
        z=float(z)/np.linalg.norm(z)

    return np.divide(v,np.asarray([x,y,z]))

def rigidTransform(A,B):

    '''map coordinates of A onto those of B'''

    assert len(A)==len(B)

    # total number of points
    N=A.shape[0];

    # find the centroids i.e. the geometric mean of each set of coordinates
    cenA=np.mean(A,axis=0)
    cenB=np.mean(B,axis=0)

    # centre the points to remove the translation component
    AA=A-np.tile(cenA,(N,1))
    BB=B-np.tile(cenB,(N,1))

    # find the covariance matrix H
    H=np.transpose(AA)*BB

    # SVD will decompose a matrix into three other matrices such that [U,S,V]=SVD(E) and E=USV**T
    U,S,Vt=np.linalg.svd(H)

    R=Vt.T*U.T

    # special reflection case
    if np.linalg.det(R)<0:
        print('Reflection detected')
        Vt[2,:]*=-1
        R=Vt.T*U.T

    # find the translation t
    t=-R*cenA.T+cenB.T

    # return what we have found
    return R,t

def findSites(prot,noncanonical=False,cull=True):

    '''locate the potential N-linked glycosylation sites (in NxS/T[/C] sequons, x!=P; excluding adjacent pairs of sequons if cull is True), returning their resnums'''

    # allow for noncanonical glycosylation sites
    if noncanonical==False:
        triplet3=['T','S']
    else:
        triplet3=['T','S','C']

    # exclude non-surface residues
    subprocess.Popen('echo q | gmx make_ndx -f Protein.pdb',shell=True).wait()
    f=open('index.ndx','r')
    indexTop=f.readlines()
    f.close()
    end=int(indexTop[-1].split()[-1])
    indexTop.append('[ resX ]\n')
    areaList=[]
    iterRes=iter(prot.residues[:-1])
    for r in iterRes:
        newLine=' '
        for a in r.atoms:
            if a.id>end:
                next(iterRes)
            else:
                newLine+='%s '%(a.id)
        newLine+='\n'
        indexTop.append(newLine)
        f=open('resi.ndx','w')
        for l in indexTop:
            f.write(l)
        f.close()
        indexTop.pop(-1)
        os.system('gmx sasa -s Protein.pdb -surface protein -output resX -n resi.ndx -probe 0.30')
        f=open('area.xvg','r')
        linesArea=f.readlines()
        f.close()
        areaList.append('%s,%s\n'%(r.resnum,linesArea[-1].split()[-1]))
        os.system('rm area.xvg')

    surfaceResi=[]
    iterArea=iter(areaList)
    for i in iterArea:
        if float(i.split(',')[1])<=0.3:
            next(iterArea)
        else:
            surfaceResi.append(int(i.split(',')[0]))

    # find the PNGs
    siteList=[]
    iterSites=iter(range(len(prot.residues.sequence())-2))
    for i in iterSites:
        if i not in surfaceResi:
            next(iterSites)
        if prot.residues.sequence()[i]=='N' and prot.residues.sequence()[i+1]!='P' and prot.residues.sequence()[i+2] in triplet3:
            siteList.append(i+prot.residues.resnums[0])

    # now we are going to remove glycosites that are immediately adjacent to each other
    if cull==True:
        for i in range(len(siteList)-1):
            if siteList[i+1]-siteList[i]==1:
                siteList.pop(i)

    return siteList

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

def isContact(r,ag,cutoff=4.0): ###will need this, but not sure if it's actually being used properly (or at all) in the current attempt... ### TOFIX

    '''determine whether or not there are direct contacts between atoms in residue r and atomgroup ag'''

    contact=False
    for a in ag.atoms:
        for rA in r.atoms:
            if np.linalg.norm(a.position-rA.position)<=cutoff:
                iscontact=True
    return contact


# set stuff up

coreDict={
        'N':['Asn-GlcNAc-4','GlcNAc-GlcNAc-4','GlcNAc-Man-4','Man-Man-3','Man-Man-6'],
        'O1':['Ser-GalNAc-4','GalNAc-Gal-3'],
        'O2':['Ser-GalNAc-4','GalNAc-Gal-3','GalNAc-GlcNAc-6'],
        'O3':['Ser-GalNAc-4','GalNAc-GlcNAc-3'],
        'O4':['Ser-GalNAc-4','GalNAc-GlcNAc-3','GalNAc-GlcNAc-6'],
        'O5':['Ser-GalNAc-4','GalNAc-GalNAc-3'],
        'O6':['Ser-GalNAc-4','GalNAc-GlcNAc-6'],
        'O7':['Ser-GalNAc-4','GalNAc-GalNAc-6'],
        'O8':['Ser-GalNAc-4','GalNAc-Gal-3'],
        'C':['Trp-Man-3']
        }
#issue with O linked and typ, therefore have a thing that will change these if the site is not Ser ### TOFIX

bondDict={
        'N':0.520,
        None:0.530 ### needs editing when more data are available ###
        }

angleDict={
        'N':120.0,
        None:110 ### Needs editing!!!
        }

mapDict={
        'GlcA':'GCA',
        'Glc':'GMY',
        'Man':'GMY',
        'Gal':'GMY',
        'GlcNAc':'GNC',
        'GlcNA':'GNC',
        'GlcNAs':'GNC',
        'GalNAc':'GNC',
        'GalNAc':'GNC',
        'Neu5Ac':'NMC',
        'Neu5Gc':'NMC',
        'Fuc':'FUC',
        'Rha':'FUC',
        'Xyl':'XYL'
        }


resNameDict={
        'Man':['xkcd:green','o'],
        'Glc':['xkcd:blue','o'],
        'Gal':['xkcd:yellow','o'],
        'GlcNAc':['xkcd:blue','s'],
        'GalNAc':['xkcd:yellow','s'],
        'GlcA':['xkcd:cyan','D'],
        'GlcN':['xkcd:cyan','o'],
        'Neu5Ac':['xkcd:purple','D'],
        'Fuc':['xkcd:red','^'],
        'Xyl':['xkcd:orange','*']
        }
############################################ dictionaries and graphs need to be looked into in considerably more detail than currently the case


# get user input

parser=argparse.ArgumentParser()
parser.add_argument('-p','--prot',default='Protein',help='prefix of the itp and pdb files for the CG protein [default: Protein]')
parser.add_argument('-c','--coreType',default='N',help='nature of the glycan core i.e. N,O,C,S,IC [default: N]')
parser.add_argument('-s','--sites',default=None,nargs='+',help='sites to be glycosylated; if not provided, they will be predicted from the sequence (providing sites is necessary for all types of glycosylation other than N-linked)')
parser.add_argument('-t','--type',default='complex',help='type of glycans to be added i.e. complex, high-mannose (only relevant for N-linked)')
parser.add_argument('-g','--glycosaminoglycan',default=None,help='glycosaminoglycan to be added [options: chondroitin sulphate, keratan sulphate, heparan sulphate, dermatan sulphate]')
parser.add_argument('-m','--monosaccharide',default=None,help='monosaccharide to be attached')
parser.add_argument('-o','--outName',default='glycoprotein',help='name for output files [default: glycoprotein]')
parser.add_argument('--bilayer',action='store_true',default=False,help='is there a bilayer present in the pdb file? [default: no]')
parser.add_argument('--viral',action='store_true',default=False,help='is the protein viral? [default: no]')
parser.add_argument('--cull',action='store_true',default=True,help='only glycosylate one site where two are adjacent to each other? [default: yes]')
parser.add_argument('--noncanonical',action='store_true',default=False,help='include noncanonical sequons? [default: no]')
args=parser.parse_args()


# get the protein coordinates, topology etc.

print('prepping to glycosylate')

u=MDAnalysis.Universe('%s.pdb'%(args.prot))
prot=u.select_atoms('protein')

if args.sites:
    sites=args.sites
else:
    sites=list(set(findSites(prot,noncanonical=args.noncanonical,cull=args.cull)))
siteList=sorted([int(x) for x in sites])
sites=siteList


# if a bilayer is present, make sure we don't build glycans into the bilayer
if args.bilayer==True:
    v=MDAnalysis.Universe('bilayer.pdb')
    bilayer=v.select_atoms('not protein and not name W and not name WN and not resname ION and not name NA+ and not name CL- and not name NA and not name CL and not name Na and not name Cl and not name Na+ and not name Cl- and not resname ION')
    tails=bilayer.select_atoms('name C* or name D*')
    for s in sites:
        r=prot.residues[s-1]
        if isContact(r,tails):
            sites.pop(sites.index(s))
       # need to account for topology somehow #####################################!!!####################################################################
    outResi=[]
    inResi=[]
    tmResi=[]
    phos=bilayer.select_atoms('name PO4')


    leaflets=LeafletFinder(v,'name PO4')
    leaflet0=leaflets.groups()[0]
    leaflet1=leaflets.groups()[1]
    if leaflet0.center_of_mass()[2]>leaflet1.center_of_mass()[2]:
        leaflet0=leaflets.groups()[1]
        leaflet1=leaflets.groups()[0]
    for r in prot.residues:
        if np.mean(r.atoms.positions)<leaflet0.center_of_mass()[2]:
            inResi.append((r.resnum,r.resname))
        if np.mean(r.atoms.positions)>leaflet1.center_of_mass()[2]:
            outResi.append((r.resnum,r.resname))
        if np.mean(r.atoms.positions)<leaflet1.center_of_mass()[2] and np.mean(r.atoms.positions)>leaflet0.center_of_mass()[2]:
            tmResi.append((r.resnum,r.resname))

    countBasicIn,countBasicOut=0,0
    for r in inResi:
        if r[1] in ['ARG','LYS','Arg','Lys']:
            countBasicIn+=1
    for r in outResi:
        if r[1] in ['ARG','LYS','Arg','Lys']:
            countBasicOut+=1
    if countBasicOut>countBasicIn:
        Tmp=inResi
        inResi=outResi
        outResi=Tmp
    for s in sites:
        for r in inResi:
            if int(s)==int(r[0]):
                sites.pop(sites.index(s))
        for r in tmResi:
            if int(s)==int(r[0]):
                sites.pop(sites.index(s))

print('glycosylation sites are: %s'%sites)

# then we run the predictor as required############################################################################


################################### NEEDS WORK ###################################
if args.type==None:
    glycanType='cpx' ### this will be a predictor applied to each site in turn, then have a list that corresponds to the sites themselves
elif args.type in ['complex','cpx']:
    glycanType='cpx'
elif args.type in ['high-mannose','hmn','oligomannose']:
    glycanType='hmn'
elif args.type in ['hybrid','hyb']:
    glycanType='hyb'
else:
    glycanType='mono'

core=args.coreType
if core=='O':
    core=random.choice(['O1','O2','O3','O4'])

# prepare the topology file #
outTop=open('topol.top','w')
outTop.write(
'''#include "martini_v3.0.3.itp"
#include "Protein.itp"
#ifdef POSRES
#include "posre.itp"
#endif
''')

for site in sites:
    outTop.write('#include "glycan-%s.itp\n'%(site))

outTop.write(
'''#include "martini_v3.0_phospholipids.itp"
#include "martini_v3.0_solvents.itp"
#include "martini_v3.0_ions.itp"

[ system ]
; name
Glycoprotein

[ molecules ]
; name  number
Protein 1
''')

for site in sites:
    outTop.write('Glycan-%s 1\n'%(site))

outTop.write('\n\n[ intermolecular_interactions ]\n[ bonds ]\n')

# Onward!

print('Doing sites:\n')
print(sites)

glycanbead=prot.residues[-1].atoms.ids[-1]+1

iterSites=iter(sites)

for site in iterSites:

    print('doing site %s'%(site))

    # lets set up for the topology
    sidebead=prot.residues[int(site)-1].atoms.ids[-1]

    # here we need to check for enough space
    ######
    # perhaps we can see if there are any other protein residues within, or check if it's solvent accessible?i
    clash=False
    iterRes=iter(prot.residues)
    for tmpRes in iterRes:
        if tmpRes.resnum in [site-1,site+1,site]:
            next(iterRes)
        else:
            separation=np.linalg.norm(cdist(tmpRes.atoms.positions,prot.residues[int(site)-1].atoms.positions))
            if separation<=10.0:
                clash=True
            else:
                clash=False
    if clash==True:
        next(iterSites)

    outTop.write('%s %s   6     0.42   100\n'%(sidebead,glycanbead))

    # before anything else, we should build the glycan
    subprocess.Popen('python build-glycans.py -t %s --%s -o glycan-%s'%(core,glycanType,site),shell=True).wait()

    ID=prot.residues[int(site)-1].atoms.ids[-1] # PFIX
    branchRes=[]
    i=1
    if len(coreDict[core])>=2:
        branchRes.append(int(site)+len(coreDict[core])-1+i)
        i-=1
    addID=prot.residues[int(site)-1].atoms.ids[-1]
    addRes=prot.residues[int(site)-1].resnum
    v=MDAnalysis.Universe('glycan-%s.pdb'%(site))
    agCore=v.select_atoms('all')
    newAtoms=[]

    siteID=ID

    pdbCore=MDAnalysis.Universe('glycan-%s.pdb'%(site))
    agCore=pdbCore.select_atoms('all')

    coords1=agCore.atoms.positions[:2]
    coords2=agCore.atoms.positions
    coords0=prot.residues[site-1].atoms.positions

    rSele=''
    resi=prot.residues[site-1]
    print(resi.resname)
    for r in prot.residues:
        d=(np.linalg.norm(r.atoms.center_of_geometry()-resi.atoms.center_of_geometry()))
        if d<=10:
            rSele+='%s '%(r.resnum)
    group=prot.select_atoms('resnum %s'%(rSele))

    v1=group.atoms.positions[0]-group.atoms.positions[1]
    v1/=unit(v1)
    newCoords=coords0+v1#(v1*1.5) ################### AHOY

    vector=group.atoms.positions[1]-group.atoms.positions[0]
    normalisedVector=vector/np.linalg.norm(vector)
    newCoords=np.asarray([newCoords[0]+(0.55*normalisedVector),newCoords[1]+(0.55*normalisedVector)])

    A=np.matrix(coords1[:2])
    B=np.matrix(newCoords[:2]) # transforming the glycan onto the side chain, based on newCoords
    R,t=rigidTransform(A,B)
    A2=(R*A.T)+np.tile(t,(1,A.shape[0]))
    A2=A2.T
    C=np.matrix(coords2)
    C2=(R*C.T)+np.tile(t,(1,C.shape[0]))
    C2=C2.T
    C2=np.asarray(C2,dtype=np.float64)
    C2=np.around(C2,decimals=3)
    agCore.atoms.positions=C2
    if np.sum(cdist(C2,prot.atoms.positions)<4.0) >=5:
        next(iterSites)

    f=open('tmp.pdb','w')

    for a in prot.atoms:
        f.write('{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'.format('ATOM',a.id,a.name,' ',a.resname,' ',a.resnum,' ',a.position[0],a.position[1],a.position[2],1.00,0.00,' ',' ')) # write everything before the addition

    for a in agCore.atoms:
        f.write('{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'.format('ATOM',a.id+ID,a.name,' ',a.resname,' ',a.resnum+addRes,' ',a.position[0],a.position[1],a.position[2],1.00,0.00,' ',' '))

    f.close()

    subprocess.Popen('gmx editconf -f tmp.pdb -o renumbered.pdb -resnr 1',shell=True).wait()
    subprocess.Popen('mv renumbered.pdb tmp.pdb',shell=True).wait()

    u=MDAnalysis.Universe('tmp.pdb')
    prot=u.select_atoms('all')
    glycanbead=prot.residues[-1].atoms.ids[-1]+1

subprocess.Popen('mv tmp.* %s.pdb'%(args.outName),shell=True).wait()
outTop.close()

print('He is all my art now\n')
