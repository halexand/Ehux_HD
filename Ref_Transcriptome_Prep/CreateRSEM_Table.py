#!/usr/bin/env python

#A script to parse out sequences longer than a certain, user defined threshhold. 

from Bio import SeqIO
import sys

def createRSEMTable(inTrans, columnName=None):
    #input the transcriptome (MMETSP) and the column within the description that you want to be the "gene" grouping. 
    inT=open(inTrans, 'rU')
    Hash={}
    if columnName:
        outFile=open(inTrans+'.RSEM.tab', 'w')

    for P in SeqIO.parse(inT, 'fasta'):
        a={}
        for i in P.description.split():
            i=i.strip()
            if '=' in i:
                j=i.split('=')
                a[j[0]]=j[1]
        Hash[P.id]=a
        if columnName:
            outFile.write(Hash[P.id][columnName])
            outFile.write('\t')
            outFile.write(P.id)
            outFile.write('\n')
    if columnName:        
        outFile.close()
    else:
        print 'Choose the variable that you would like to be used for RSEM gene families. Choices are as follows:'
        for l in Hash[P.id].keys():
            print l


if __name__=="__main__":
    print 'Loading transcriptomes'
    createRSEMTable('Ehux_All_Transcripts_Cleaned_RSEM.nt.fa', '/OrthGrp')
    
    
