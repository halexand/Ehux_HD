#/usr/bin/env python
#Read in important paths 
import pandas as pd
import cPickle as cpk
import matplotlib as mpl
import sys
import numpy as np
sys.path.append("/Users/harrietalexander/anaconda/lib/python2.7/site-packages/matplotlib_venn-0.11-py2.7.egg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from matplotlib import gridspec
from itertools import combinations
import palettable.colorbrewer as b2m
import palettable as pal
import glob
from Bio import SeqIO
import os


def Calculate_Venn_comparison(hash):
    #input = a Hash of the form hash[Organism]=list of orthologus groups
    #Outputs a hash of the form key = all the variations of thes sets : list of orthologus genes in an inclusive way
    #It is inclusive so the set containing A uniq genes is going to contain all the genes for an org
    variations={}
    for i in range(len(hash)):
        for v in combinations(hash.keys(), i+1):
            vsets = [hash[x] for x in v]
            variations[tuple(sorted(v))]=reduce(lambda x,y: x.intersection(y), vsets)
    return variations

def GetDifference_Venn(variations):
    #Function takes the input from the above function and outputs a hash of the same format but 
    #each set is uniqe: e.g. there are no repeats of orthologus groups across sets. 
    outdict={}
    vkeys=variations.keys()
    #loop over each of the variable cases
    for v in variations:
        #create a set to do the comparisons
        vset=set(v)
        vdata_set=set(variations[v])
        #loop over all other variations
        for j in variations:
            #if v is a subset of jset we want to remove the items of jset from v set
            jdata_set=set(variations[j])
            jset=set(j)
            if vset.issubset(jset):
                if vset==jset:
                    pass
                else: 
                    newdata=vdata_set-jdata_set
                    vdata_set=newdata
        outdict[v]=vdata_set
    return outdict

def loadData(inFile):
    #load data from orthoMCL
    for i,c in enumerate(inFile):
#     if i==10:
#         break
        if i%10000 == 0:
            print i
        c=c.split()
        geneFamily=c[0].strip()
        transcript=c[1].strip()
        if geneFamily in gFHash.keys():
            gFHash[geneFamily].append(transcript)
        else:
            gFHash[geneFamily]=[transcript]
    cpk.dump(gFHash, open('Ehux_Dictionary.pickle', 'w'))
    return gFHash


if __name__=='__main__': 
    if os.path.isfile('Ehux_Dictionary.pickle'):
        gFHash=cpk.load(open('Ehux_Dictionary.pickle', 'r'))
    else:
        inFile=open("Ehux_All_Transcripts_Cleaned_RSEM.nt.fa.RSEM.tab", 'r')
        loadData(inFile)
#Parse the gFHash to tally the number of genes in an orthologus group and the number of genes from each of the individual taxa
#Panda dataframe of form : Orthologus group | total number of genes in orthologus group | total number from each strain
#Nested for loop... so it takes a while. 
    
    Hist_PD=pd.DataFrame(index=gFHash.keys(),columns=['NumGenes', 'Emi374', 'Emi379', 'Emi370', 'Emi219', 'Emihu1'])
    Hist_PD=Hist_PD.fillna(0)
    for i,key in enumerate(gFHash):
        l=len(gFHash[key])
        Hist_PD.loc[key, 'NumGenes']=l
        for transcript in gFHash[key]:
            org=transcript.split('|')[0]
            Hist_PD.loc[key,org]+=1
    GenesInOrg={}
    GenesInOrg['Emi219']=Hist_PD[Hist_PD.Emi219>0].index
    GenesInOrg['Emi379']=Hist_PD[Hist_PD.Emi379>0].index
    GenesInOrg['Emi370']=Hist_PD[Hist_PD.Emi370>0].index
    GenesInOrg['Emi374']=Hist_PD[Hist_PD.Emi374>0].index
    GenesInOrg['Emihu1']=Hist_PD[Hist_PD.Emihu1>0].index

    Genes_In_Each_Cat=Calculate_Venn_comparison(GenesInOrg)
    Genes_Uniq_Venn=GetDifference_Venn(Genes_In_Each_Cat)
#write out lists of each of the different gene sets.
    print Genes_Uniq_Venn
    for key in Genes_Uniq_Venn:
        fName='_'.join(list(key))
        fName='Species_GeneLists_Venn/'+fName+'.list'
        file=open(fName, 'w')
        for item in Genes_Uniq_Venn[key]:
            file.write(item.strip())
            file.write('\n')
        file.close()
