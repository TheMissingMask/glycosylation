#!/usr/bin/python

import argparse
import MDAnalysis
import glob
import os

aaDict={
        'ARG':'R',
        'HIS':'H',
        'LYS':'K',
        'ASP':'D',
        'GLU':'E',
        'SER':'S',
        'THR':'T',
        'ASN':'N',
        'GLN':'Q',
        'CYS':'C',
        'GLY':'G',
        'PRO':'P',
        'ALA':'A',
        'ILE':'I',
        'LEU':'L',
        'MET':'M',
        'PHE':'F',
        'TRP':'W',
        'TYR':'Y',
        'VAL':'V'
        }

mapDict={
        'O':0,
        'i':1,
        'M':2,
        'o':3
        }

parser=argparse.ArgumentParser()
parser.add_argument('-p','--prot')
args=parser.parse_args()

u=MDAnalysis.Universe(args.prot)
prot=u.select_atoms('protein')

chainIDs=list(set(prot.chainIDs))
if len(chainIDs)>0:
    for chain in chainIDs:
        prot=u.select_atoms('protein and segid %s'%(chain))

        seq=[]
        fasta='>protein sequence\n'
        for r in prot.residues.resnames:
            fasta+=aaDict[r]
            seq.append(aaDict[r])

        f=open('prot-%s.fasta'%(chain),'w')
        for l in fasta:
            f.write(l)
        f.close()

        os.system('tmhmm -f prot-%s.fasta -m TMHMM2.0.model -p'%(chain))

        topology=[]
        txt=glob.glob('protein.annotation')
        f=open(txt[0])
        lines=f.readlines()
        f.close()
        for l in lines[1:]:
            for k in l.strip():
                topology.append(k)

        ca=prot.select_atoms('name CA')
        resnums=ca.residues.resnums

        bfactors=[]
        for i in range(len(resnums)):
            bfactors.append(mapDict[topology[i]])

        rbfactors=[]
        for i in range(len(prot.residues)):
            b=bfactors[i]
            n=len(prot.residues[i].atoms)
            for j in range(n):
                rbfactors.append(b)
        prot.atoms.tempfactors=rbfactors
        prot.write('bfactors-prot%s.pdb'%(chain))

else:

    seq=[]
    fasta='>protein sequence\n'
    for r in prot.residues.resnames:
        fasta+=aaDict[r]
        seq.append(aaDict[r])

    f=open('prot.fasta','w')
    for l in fasta:
        f.write(l)
    f.close()

    os.system('tmhmm -f prot.fasta -m TMHMM2.0.model -p')

    topology=[]
    txt=glob.glob('protein.annotation')
    f=open(txt[0])
    lines=f.readlines()
    f.close()
    for l in lines[1:]:
        for k in l.strip():
            topology.append(k)

    ca=prot.select_atoms('name CA')
    resnums=ca.residues.resnums

    bfactors=[]
    for i in range(len(resnums)):
        bfactors.append(mapDict[topology[i]])

    rbfactors=[]
    for i in range(len(prot.residues)):
        b=bfactors[i]
        n=len(prot.residues[i].atoms)
        for j in range(n):
            rbfactors.append(b)
    prot.atoms.tempfactors=rbfactors
    prot.write('bfactors-prot.pdb')
