#!/usr/bin/env python

#A script to parse out sequences longer than a certain, user defined threshhold. 

from Bio import SeqIO
import sys

def prepTranscriptomes(inTrans, inProteome, inPassProtein):
    #input the transcriptome (MMETSP), the proteom (MMETSP), and the orthoMCL compliantFasta file 
    inT=open(inTrans, 'rU')
    inPP=open(inPassProtein, 'rU')
    inP=open(inProteome, 'rU')
    P_pass=[]
    for P in SeqIO.parse(inPP, 'fasta'):
        d=P.id.split('|')#get the passing gene ids from proteins
        P_pass.append(d[1])
        specName=d[0]
    proteinHash={}    

    for P in SeqIO.parse(inP, 'fasta'): #get the nuc name that corresponds to the protein name
        d=P.description.split("/")
        proteinHash[specName+'|'+d[0].strip()]=d[5].strip().split('=')[1]
    
#    PassingNucleotides=set()
#    for p in P_pass:
#        n=[proteinHash[p].split("=")[1]]
#        PassingNucleotides.update(n)
    return proteinHash

def importGenome(genome_transcripts):
    iG=open(genome_transcripts, 'rU')
    outHash={}
    for i in SeqIO.parse(iG, 'fasta'):
        name=i.id
        sname="|".join(name.split('|')[1:3])
        outHash[sname]=name
    return outHash

def importOrthMCL(groups, singletons):
    #create a hash containgin the groups and singletons from orthoMCL clustering; we will be using this to add group names etc. to the names and ultimately to the description/name of the fasta used for mapping
    iG=open(groups, 'rU')
    iS=open(singletons,'rU')
    orthMCLgroups={}
    num_=0
    for g in iG:
        gs=g.strip().split()
        groupName=gs[0].strip(":")
        num=int(groupName.split("_")[2])
        if num>num_:
            num_=num
        orthMCLgroups[groupName]=gs[1:]
   
    preName='_'.join(groupName.split("_")[:2])
    for s in iS:
        s=s.strip()
        num_+=1
        newGrp=preName+"_"+str(num_)
        orthMCLgroups[newGrp]=[s]
    return orthMCLgroups

#get the proteins and assign "gene familiy" names to them

def invert(d):
    return dict( (v,k) for k in d for v in d[k] )
def invertDictionary(orig_dict):
    result = {} # or change to defaultdict(list)
    for k, v in orig_dict.iteritems():
        result.setdefault(v, []).append(k)
    return result
def merge_dicts(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    '''
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


if __name__=="__main__":
    print 'Loading transcriptomes'
    Ehu1=importGenome('Emihu1_reduced_transcripts.fasta')
    E374=prepTranscriptomes('Emiliania-huxleyi-374.nt.fa', 'Emiliania-huxleyi-374.pep.fa', 'Emi374.fasta')
    E379=prepTranscriptomes('Emiliania-huxleyi-379.nt.fa', 'Emiliania-huxleyi-379.pep.fa', 'Emi379.fasta')
    E370=prepTranscriptomes('Emiliania-huxleyi-CCMP370.nt.fa', 'Emiliania-huxleyi-CCMP370.pep.fa', 'Emi370.fasta')
    E219=prepTranscriptomes('Emiliania-huxleyi-PLYM219.nt.fa', 'Emiliania-huxleyi-PLYM219.pep.fa', 'Emi219.fasta')
    
    orthoMCLGroups=importOrthMCL('named_groups_1.5.txt', 'named_singletons.txt')
    orthoMCLGroups_inverted=invert(orthoMCLGroups)
    AllDict=merge_dicts(Ehu1, E374, E379, E370, E219)
    AllDict_inverted=invertDictionary(AllDict)
    newFasta=open('output_fasta.fa', 'w')
    ortho_inverted_set=set(orthoMCLGroups_inverted.keys())
    outSeq=[]
    c=0
    print 'Preparing output fasta file'
    for seq in SeqIO.parse(open('AllEhu_Transcripts.nt.fa', 'rU'), 'fasta'):
        c+=1
        n=seq.id
        if n in AllDict_inverted.keys():
            k=AllDict_inverted[n][0]
            if k in orthoMCLGroups_inverted.keys():   
                seq.description=seq.description+" /OrthGrp="+ orthoMCLGroups_inverted[k]+" /ProteinID="+k
                seq.id=k
                outSeq.append(seq)
        if c%1000==0:
            print c
    SeqIO.write(outSeq, newFasta, 'fasta')
    
