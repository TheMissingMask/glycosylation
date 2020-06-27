#!/usr/bin/python

import os
import MDAnalysis
import argparse


parser=argparse.ArgumentParser()
parser.add_argument('-p','--prot')
parser.add_argument('-r','--radius',default='0.782')
args=parser.parse_args()

os.system('echo "1\n" | gmx sasa -oa -or -o -probe %s -s %s'%(args.radius,args.prot))

cutoff=float(args.radius)

f=open('atomarea.xvg')
lines=f.readlines()
f.close()
bfactors=[]
for l in lines:
    cols=l.split()
    if cols[0][0] in ['@','#']:
        next
    else:
        if float(cols[1])>=cutoff:
            bfactors.append(1)
        else:
            bfactors.append(0)
            
u=MDAnalysis.Universe(args.prot)
prot=u.select_atoms('protein')
ca=prot.select_atoms('name CA')
bfactorsRes=bfactors

prot.tempfactors=bfactorsRes
prot.write('bfactor-sasa-%s.pdb'%(args.radius))

f=open('bfactor-sasa-%s.pdb'%(args.radius))
lines=f.readlines()
f.close()

buriedLines=[]
surfaceLines=[]
for l in lines:
    cols=l.split()
    if cols[0]!='ATOM':
        next
    else:
        b=cols[-3]
        if float(b)==0.00:
            buriedLines.append(l)
        elif float(b)==1.00:
            surfaceLines.append(l)

f=open('surface.pdb','w')
for l in surfaceLines:
    f.write(l)
f.close()
f=open('buried.pdb','w')
for l in buriedLines:
    f.write(l)
f.close()
