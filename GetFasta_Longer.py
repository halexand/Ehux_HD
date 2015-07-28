#!/usr/bin/env python

#A script to parse out sequences longer than a certain, user defined threshhold. 

from Bio import SeqIO
import sys

def getSeqLonger(inFile, cutoff):
    inF=open(inFile, 'rU')
    LongSeq=[]
    ShortSeq=[]
    for record in SeqIO.parse(inF, "fasta"):
        if len(record.seq)>int(cutoff):
            LongSeq.append(record)
        else: 
            ShortSeq.append(record)
    Sout=open(inFile+".ShorterThan"+str(cutoff),'w')
    Lout=open(inFile+".LongerThan"+str(cutoff),'w')
    SeqIO.write(ShortSeq, Sout, 'fasta')
    SeqIO.write(LongSeq, Lout, 'fasta')
    Lout.close()
    Sout.close()
    print "greater than ", cutoff, ": ", str(len(LongSeq))     
    print "less than ", cutoff, ": ", str(len(ShortSeq))     

if __name__=="__main__":
    getSeqLonger(sys.argv[1], sys.argv[2])
