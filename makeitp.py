#!/usr/bin/python

import argparse
import glob

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

parser=argparse.ArgumentParser()
parser.add_argument('-p','--prot',default='Protein',help='prefix of the itp and pdb files for the CG protein [default: Protein]')
parser.add_argument('--sites',nargs='+')
args=parser.parse_args()



glycanList=[]
#sites=[]
glycanFiles=glob.glob('glycan-*.itp')
for f in glycanFiles:
    glycanList.append(f.split('.')[0])
#sites=args.sites
#glycanList=['glycan-%s'%(x) for x in sites]
print(glycanList)

#itpLines=linesFromFile('%s.itp'%(args.prot))
itpLines=linesFromFile('Protein.itp')


for glycan in glycanList:

    maxID=int(linesBetween(itpLines,'[ atoms ]','[ bonds ]')[-1].split()[0])
    maxRes=int(linesBetween(itpLines,'[ atoms ]','[ bonds ]')[-1].split()[2])

    itpDict={
    '[ atoms ]':[],
    '[ bonds ]':[],
    '[ constraints ]':[],
    '[ angles ]':[],
    '[ dihedrals ]':[],
    '[ exclusions ]':[],
    '[ virtual_sitesn ]':[]
    }
    glycanLines=linesFromFile('%s.itp'%(glycan))

    idMap={}
    resMap={}

    for l in linesBetween(glycanLines,'[ atoms ]','[ bonds ]'):
        cols=l.split()
        atomid=int(cols[0])
        resid=int(cols[2])
        idMap[atomid]=atomid+maxID
        resMap[resid]=resid+maxRes
    
    for l in linesBetween(itpLines,'[ atoms ]','[ bonds ]'):
        itpDict['[ atoms ]'].append(l)
    
    for l in linesBetween(glycanLines,'[ atoms ]','[ bonds ]'):
        cols=l.split()
        newLine='%s  %s  %s  %s  %s   %s  %s\n'%(idMap[int(cols[0])],cols[1],resMap[int(cols[2])],cols[3],cols[4],idMap[int(cols[5])],cols[6])
        itpDict['[ atoms ]'].append(newLine)
    
    for l in linesBetween(itpLines,'[ bonds ]','[ constraints ]'):
        itpDict['[ bonds ]'].append(l)
    
    for l in linesBetween(glycanLines,'[ bonds ]','[ constraints ]'):
        cols=l.split()
        newLine='%s  %s  %s  %s  %s\n'%(idMap[int(cols[0])],idMap[int(cols[1])],cols[2],cols[3],cols[4])
        itpDict['[ bonds ]'].append(newLine)
    
    addResid='%s'%(glycan).split('-')[-1]
    for l in linesBetween(itpLines,'[ atoms ]','[ bonds ]'):
        cols=l.split()
        if int(cols[2])==int(addResid) and cols[4]=='SC1':
            bondID=int(cols[0])
        
    newLine='%s  %s  %s  %s  %s\n'%(idMap[1],bondID,2,0.40,10)
    itpDict['[ bonds ]'].append(newLine)

    for l in linesBetween(itpLines,'[ constraints ]','[ angles ]'):
        itpDict['[ constraints ]'].append(l)
    
    for l in linesBetween(glycanLines,'[ constraints ]','[ angles ]'):
        cols=l.split()
        newLine='%s  %s  %s  %s\n'%(idMap[int(cols[0])],idMap[int(cols[1])],cols[2],cols[3])
        itpDict['[ constraints ]'].append(newLine)
    
    for l in linesBetween(itpLines,'[ angles ]','[ dihedrals ]'):
        itpDict['[ angles ]'].append(l)
    
    for l in linesBetween(glycanLines,'[ angles ]','[ dihedrals ]'):
        cols=l.split()
        newLine='%s  %s  %s %s %s %s\n'%(idMap[int(cols[0])],idMap[int(cols[1])],idMap[int(cols[2])],cols[3],cols[4],cols[5])
        itpDict['[ angles ]'].append(newLine)
    
    for l in linesBetween(itpLines,'[ dihedrals ]','[ exclusions ]'):
        itpDict['[ dihedrals ]'].append(l)
    
    for l in linesBetween(itpLines,'[ virtual_sitesn ]',''):
        itpDict['[ virtual_sitesn ]'].append(l)
    
    for l in linesBetween(glycanLines,'[ virtual_sitesn ]',''):
        cols=l.split()
        print(l)
        newLine='%s  %s  %s %s %s\n'%(idMap[int(cols[0])],cols[1],idMap[int(cols[2])],idMap[int(cols[3])],cols[4])
        itpDict['[ virtual_sitesn ]'].append(newLine)
    
    f=open('glycoprotein.itp','w')
    f.write('[ moleculetype ]\nGlycoprotein   3\n\n')
    for k in ['[ atoms ]','[ bonds ]','[ constraints ]','[ angles ]','[ dihedrals ]','[ exclusions ]','[ virtual_sitesn ]']:
        f.write('%s\n'%(k))
        for l in itpDict[k]:
            f.write(l)
        f.write('\n')
    f.close()

    itpLines=linesFromFile('glycoprotein.itp')

print('Perhaps this has worked...')
