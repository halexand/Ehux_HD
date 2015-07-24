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


def readInSpeciesVenn(directory):
    #input file to be read in    #Outputs a hash of the form key = all the variations of thes sets : list of orthologus genes in an inclusive way
    #It is inclusive so the set containing A uniq genes is going to contain all the genes for an org
    outHash={}
    
    for fin in glob.glob(directory+'/*'):
        name=fin.split('/')
        name=name[1].split('.')[0]
        a=set()
        for line in open(fin):
            line=line.strip()
            a.add(line)
        outHash[name]=a
    return outHash

def readInTranscriptGene(file):
    #read in the transcript :: OG tab file
    #create hash of form OG :: transcript
    outHash={}
    for line in open(file,'r'):
        a=line.split()
        a1=a[1].strip()
        a2=a[0].strip()
        if a1 in outHash.keys():
            outHash[a1].append(a2)
        else:
            outHash[a1]=a2
    return outHash


def readInKEGG(inFile):
    #read in KAAS annotation transcript :: Hash
    outHash={}
    for line in open(inFile,'r'):
        a=line.split()
        if len(a)>1:
            if a[0].strip() in outHash.keys():
                outHash[a[0].strip()].append(a[1].strip())       
            else:
                outHash[a[0].strip()]=a[1].strip()
    return outHash

def CreateKEGG_Hash_OG(KAAS_Hash,gFHash):
    newOG_KEGGHash={}
    for gene in gFHash:
        newOG_KEGGHash[gene]=[]
        for transcript in gFHash[gene]:
            kegg=KAAS_Hash[transcript]
            if kegg in newOG_KEGGHash[gene]:
                pass
            elif kegg==None:
                pass
            else: 
                newOG_KEGGHash[gene].append(kegg)#
    return newOG_KEGGHash


if __name__=='__main__': 
    #read in all the files
    testDir='test'
    Transcript_To_OG=cpk.load(open('../orthoMCL_output/Ehux_Dictionary.pickle', 'r'))
    print 'reading in files from '+testDir
    Venn_Set_Hash=readInSpeciesVenn(testDir)
    #Transcript_To_OG=readInTranscriptGene('Ehux_All_Transcripts_Cleaned_RSEM.nt.fa.RSEM.tab')
    KeggHash=readInKEGG('KAAS_Ehux_Protein_Annotation.tab')
    #create the OG :: Kegg Hash
    print 'Getting KEGG Dictionary'
    Kegg_OG_Hash=CreateKEGG_Hash_OG(KeggHash, Transcript_To_OG)
    print Kegg_OG_Hash
