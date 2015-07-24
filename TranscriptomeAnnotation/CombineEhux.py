#/usr/bin/env python
import re
import cPickle as cpk
def readInKEGG(inFile):
    #read in KAAS annotation transcript :: Hash
    outHash={}
    for line in open(inFile,'r'):
        line=line.split()
        if len(line)>1:
            outHash[line[0].strip()]=[line[1].strip()]
            
        else:
            outHash[line[0].strip()]=[]
    return outHash

def readInEhuxKegg(inFile):
    outHash={}
    for line in open(inFile, 'r'):
        if line.startswith('D'):
            l=line.strip().split()
            key=l[1]
            s=[s for s in l if re.search('K[0-9]', s)]
            ks=[k for k in s if len(k)==6]
            outHash[key]=ks
    return outHash

if __name__=='__main__': 
    KeggHash=readInKEGG('KAAS_Ehux_Protein_Annotation.tab')
    EhxKegg=readInEhuxKegg('ehx00001.keg.renamed')
    newHash=KeggHash
    for key in KeggHash:
        if key in EhxKegg.keys():
            newK=EhxKegg[key]
            if newK in KeggHash[key]:
                pass
            else:
                KeggHash[key].append(newK)
    cpk.dump(KeggHash, open('Kegg_Hash.pickle', 'w'))  

