#!/usr/bin/python

import MDAnalysis
#gmx sasa -oa -or -o -s 1mqm.pdb 

cutoff=0.025 # should it be 0.025?

#f=open('resarea.xvg')
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
            
u=MDAnalysis.Universe('1mqm.pdb')
prot=u.select_atoms('protein')
ca=prot.select_atoms('name CA')
#ca.tempfactors=bfactors
bfactorsRes=bfactors
        
#for i in range(len(bfactors)):
#    n=len(ca.residues[i].atoms)
#    b=bfactors[i]
#    for j in range(n):
#        bfactorsRes.append(b)
        
prot.tempfactors=bfactorsRes
prot.write('bfactor-sasa.pdb')
