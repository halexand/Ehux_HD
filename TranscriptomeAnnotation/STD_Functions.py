#!/usr/bin/env python

'''
Created on November 9, 2013

STD_Function.py

A generalized collection of 

@author: harrietalexander
'''

import sys
import glob
import os
import numpy as np
import re
import csv
import matplotlib.pylab as plt
import pickle

#General functions
def Write_Dict_To_File(dict,fileName,header=None):
	#Writes dictionary to file specieified header row and dictionary
	writer=csv.writer(open(fileName, 'wb'), delimiter='\t')
	writer.writerow(header)
	for name, value in dict.items():
		test = [name] + value
		writer.writerow(test, )

def makeLists(hash,N):
# Given Hash and N (range or size of data) will make a new list of list from the hash by column/list entry
# Use for plotting etc. 
	keys=[]
	newList=[[] for x in xrange(N)]
	for key in hash:
		keys.append(key)
		for x in range(0,N):
			newList[x].append(hash[key][x])
	return newList , keys
	
def flatList(L):
#Given a list of sublits this will return a "flat" list (i.e. not a list of lists)
	outList=[item for sublist in L for item in sublist]
	return outList

# Functions to upload tab delimited raw counts ; 

def importDict(inFile):
#Read in tab delimited file in the format Header// seqName / count / count / count; returns dictionary
	hash={}
	handle=open(inFile)
	reader=csv.reader(handle, delimiter='\t')
	header=next(reader, None)
	for line in reader: 
		Seq=line[0]
		Nums=[int (x) for x in line[1:]]
		hash[Seq]=Nums
	head_all=[hash,header]
	return head_all
	
def TPMCalc(inHash):
#Given a hash(seqname)=[num across conditions] will calculate the TPM (dividing each entry by the sum of the columns)	
	c=0
	newHash={}
	for key in inHash:
		if c==0:
			listSum=inHash[key]
			c+=1
		else: 
			tmp=[x+y for x,y in zip(listSum, inHash[key])]
			listSum=tmp
	for key in inHash: 
		tpm=[(float(a)/float(b))*1e6 for a,b in zip(inHash[key], listSum)]
		newHash[key]=tpm
	return newHash
	
def cutDataDict(inHash, cutoff):
#Given a cutoff value (tpm etc.) function will return only all entries for which any of the treatments are greater than the cutoff
	passHash={}
	failHash={}
	passKey=[]
	failKey=[]
	for key in inHash:
		if np.max(inHash[key])>cutoff:
			passHash[key]=inHash[key]
			passKey.append(key)
		else:
			failHash[key]=inHash[key]
			failKey.append(key)
	return passHash, passKey, failHash, failKey

def SGNCCalc(StableGenes, CountsHash):	
#Normalizes the Counts to Identified Stable genes counts across the treatments
	listCounts=[]
	for genes in StableGenes:
		listCounts.append(CountsHash[genes])
	mean_listCounts=[]
	stdev_listCounts=[]
	for x in range(len(CountsHash[genes])):
		newList=[]
		for y in range(len(listCounts)):
			newList.append(listCounts[y][x])
		mean_listCounts.append(np.mean(newList))
		stdev_listCounts.append(np.std(newList))
	newHash={}
	for key in CountsHash: 
		SGNC=[(float(a)/float(b)) for a,b in zip(CountsHash[key], mean_listCounts)]
		newHash[key]=SGNC
	return mean_listCounts, stdev_listCounts, newHash

def STD_Calc(Dif_Genes, SGNC_Hash, CNplus, CNminus, CPplus, CPminus, Cinsitu):
#Import differentially regulated genes and SGNC_hash coordinate of plus N/P and minus N/P in list, range of insitu samples; calculate the STD score based for either N or P
	STDHashP={}
	STDHashN={}
	for gene in Dif_Genes:
		insitu=[]
		for i in Cinsitu: 
			insitu.append(SGNC_Hash[gene][i])
		Np=SGNC_Hash[gene][CNplus]
		Nm=SGNC_Hash[gene][CNminus]
		N=[Np,Nm]
		Pp=SGNC_Hash[gene][CPplus]
		Pm=SGNC_Hash[gene][CPminus]
		P=[Pp,Pm]
		STDP=[]
		STDN=[]
# ## option 1: Difference of the max and min
# 		for q in insitu:
# 			if (np.max(P)-np.min(P))!=0:
# 				stdP=((q-np.min(P))/(np.max(P)-np.min(P)))
# 			else:
# 				stdP=0
# 			STDP.append(stdP)
# 			if (np.max(N)-np.min(N))!=0:
# 				stdN=((q-np.min(N))/(np.max(N)-np.min(N)))
# 			else: 
# 				stdN=0
# 			STDN.append(stdN)

## option 2: Nminus - Nplus (limited minus the addition)
		for q in insitu:
			if (Pm-Pp)!=0:
				stdP=(q-Pp)/(Pm-Pp)
			else: 
				stdP=0
			STDP.append(stdP)
			if (Nm-Np)!=0:
				stdN=(q-Np)/(Nm-Np)
			else: 
				stdN=0
			STDN.append(stdN)
		STDHashP[gene]=STDP
		STDHashN[gene]=STDN
	return STDHashP, STDHashN
	
				
#Functions to work with ASC output data 
	
def testPPnotNAN(hashseq):
	"determine if gene is greater than a specific value"
	SeqPass=[]
	for seq in hashseq:
		values=[]
		for item in hashseq[seq]:
			values.append(item[1][1:])
		if np.any(np.isnan(values)):
			cc=0
		else:
			SeqPass.append(seq)
			
	return SeqPass
	
def testPP(hashseq,cutoff):
	"determine if gene is greater than a specific value"
	SeqPass=[]
	for seq in hashseq:
		values=[]
		for item in hashseq[seq]:
			values.append(item[1][1:])
		if (np.max(values)< cutoff):
			SeqPass.append(seq)
	return SeqPass

def intersect(a,b):
	return list(set(a) & set(b))
	
def importASC(directory, tail):
#Input directory where all the txt files; input directory and tail to find 
#Splits into Incubation only and field sample sites
	os.chdir(directory)
	incuHashPP={}
	fieldHashPP={}
	incucount=0
	fieldcount=0
	Ecount=0
	EHashPP={}
	for file in glob.glob(tail):
		fileName=file[:-8]
		if re.search(r'[0-9]',fileName):
			handle=open(file, "rU")
			reader=csv.reader(handle, delimiter='\t',)
			next(reader,None)
			if fieldcount==0:
				for line in reader:
					Seq=line[0]
					lst=[float(i) for i in line[1:]] 
					test=[fileName,lst]
					fieldHashPP[Seq]=[test]
				fieldcount+=1
			else:
				for line in reader:
					Seq=line[0]
					lst=[float(i) for i in line[1:]] 
					test=[fileName,lst]
					fieldHashPP[Seq].append(test)
		else:
			handle=open(file, "rU")
			reader=csv.reader(handle, delimiter='\t',)
			next(reader,None)
			if incucount==0:
				for line in reader:
					Seq=line[0]
					lst=[float(i) for i in line[1:]] 
					test=[fileName,lst]
					incuHashPP[Seq]=[test]
				incucount+=1
			else:
				for line in reader:
					Seq=line[0]
					lst=[float(i) for i in line[1:]] 
					test=[fileName,lst]
					incuHashPP[Seq].append(test)
	return incuHashPP, fieldHashPP 
	
def subHash_fromKeyList(Hash, list):
	outHash={}
	for key in list: 
		outHash[key]=Hash[key]
	return outHash

#Differentially regulated functions
def makeDict(InFile):
	hash={}
	handle=open(InFile)
	reader=csv.reader(handle, delimiter='\t')
	next(reader, None)
	for line in reader: 
		UD=[float('nan'),float('nan')]
		Seq=line[0]
		UD[0]=float(line[2])
		UD[1]=float(line[3])
		hash[Seq]=UD
	return hash

def difRegInt(hashSeq, passSeq, cutoff):
	count=0
	Up=[]
	Down=[]
	for key in hashSeq:
		if hashSeq[key][0]>cutoff:
			Up.append(key)
		if hashSeq[key][1]>cutoff:
			Down.append(key)
	intUp=intersect(Up,passSeq)
	intDown=intersect(Down,passSeq)
	return(intUp, intDown)	
		
def difRegGenes(Directory, passGenes, cutoff):
#Identify differentially regulated genes for N and P
	os.chdir(Directory)
	#Look for differentially regulated genes
	PlusNMinusN="PPSASB_2.txt" #Change of N
	PlusPMinusP="PPSCSD_2.txt" #Change of P
	MinusNPlusP="PPSBSC_2.txt" #Should have little to no change-- as both have the addition of P
	PlusNMinusP="PPSASD_2.txt" #Should have little to no change-- as both ahve the addition of N
	MinusNMinusP="PPSBSD_2.txt" #Should have some change? 
	
	pNmN=makeDict(PlusNMinusN)
	pPmP=makeDict(PlusPMinusP)
	mNpP=makeDict(MinusNPlusP)
	pNmP=makeDict(PlusNMinusP)
	mNmP=makeDict(MinusNMinusP)
	
	
	
	pNmN_Up, pNmN_Down=difRegInt(pNmN, passGenes,cutoff)
	
	print "pNmN",len(pNmN_Up),len(pNmN_Down)
	
	pPmP_Up, pPmP_Down=difRegInt(pPmP, passGenes, cutoff)
	
	print "pPmP",len(pPmP_Up),len(pPmP_Down)
	
	mNpP_Up, mNpP_Down=difRegInt(mNpP, passGenes, cutoff)
	
	print "mNpP",len(mNpP_Up),len(mNpP_Down)
	
	pNmP_Up, pNmP_Down=difRegInt(pNmP, passGenes,cutoff)
	
	print "pNmP",len(pNmP_Up),len(pNmP_Down)
	
	mNmP_Up, mNmP_Down=difRegInt(mNmP, passGenes, cutoff)
	 
	print "mNmP",len(mNmP_Up),len(mNmP_Down)

	Ngenes=[pNmN_Up, pNmN_Down]
	Pgenes=[pPmP_Up, pPmP_Down]
	return Ngenes, Pgenes
	
	
	
def count_Quadrants(STD_Counts): 
#	Input STD Counts list (broken down by Pup, Nup, etc, etc. 
#	STD_all=[PupP_STD, PupN_STD, PdnP_STD, PdnN_STD, NupP_STD, NupN_STD, NdnP_STD, NdnN_STD]
	#Compare by Pup/Pdn etc. etc. 
	Pup=STD_Counts[0:2]
	Pdn=STD_Counts[2:4]
	Nup=STD_Counts[4:6]
	Ndn=STD_Counts[6:8]
	All_Set=[Pup, Pdn, Nup, Ndn]
	
	quad1=[[],[],[],[],[]] # N and P limited (Co-limited) ; p>1 and n>1
	quad2=[[],[],[],[],[]] # P limited; p>1 and n<0
	quad3=[[],[],[],[],[]] # N and P replete; p<0 and n<0
	quad4=[[],[],[],[],[]] # N limited 
	quadPlim=[[],[],[],[],[]] # P limited 
	quadPrep=[[],[],[],[],[]] # P repelete
	quadNlim=[[],[],[],[],[]] # 
	quadNrep=[[],[],[],[],[]]
	quadnull=[[],[],[],[],[]]
	
	count=0
	for set in All_Set:
		for key in set[0]:
			count+=1
			Pnums=set[0][key]
			Nnums=set[1][key]
			for x in np.arange(len(Pnums)):
				p=Pnums[x]
				n=Nnums[x]
				if (p>=1) & (n>=1): 
					#Quadrant1 Co-limited
					quad1[x].append(key)
				elif (p>=1) & (n<=0):
					#Quadrant2 P-limited 
					quad2[x].append(key)
				elif (p<=0) & (n<=0):
					#Quadrant3 Replete
					quad3[x].append(key)
				elif (p<=0) & (n>=1):
					#Quadrant4 N-limited
					quad4[x].append(key)
				elif (p>=1) & (n>0) & (n<1): 
					quadPlim[x].append(key)
				elif (p<=0) & (n>0) & (n<1):
					quadPrep[x].append(key)
				elif (n>=1) & (p>0) & (p<1):
					quadNlim[x].append(key)
				elif (n<=0) & (p>0) & (p<1):
					quadNrep[x].append(key)
				else: 
					quadnull[x].append(key)
	allquad=[[quad1, quad2, quad3, quad4, quadPlim, quadPrep, quadNlim, quadNrep, quadnull], count]
	print np.sum([len(quad1[0]),len(quad2[0]), len(quad3[0]), len(quad4[0]), len(quadPlim[0]), len(quadPrep[0]), len(quadNlim[0]), len(quadNrep[0]), len(quadnull[0])])							
	print "quad1:", len(quad1[0]), len(quad1[1]), len(quad1[2]), len(quad1[3]), len(quad1[4]) 
	print "quad2:", len(quad2[0]), len(quad2[1]), len(quad2[2]), len(quad2[3]), len(quad2[4]) 
	print "quad3:", len(quad3[0]), len(quad3[1]), len(quad3[2]), len(quad3[3]), len(quad3[4]) 
	print "quad4:", len(quad4[0]), len(quad4[1]), len(quad4[2]), len(quad4[3]), len(quad4[4]) 
	print "quadPlim:", len(quadPlim[0]), len(quadPlim[1]), len(quadPlim[2]), len(quadPlim[3]), len(quadPlim[4]) 
	print "quadPrep:", len(quadPrep[0]), len(quadPrep[1]), len(quadPrep[2]), len(quadPrep[3]), len(quadPrep[4]) 
	print "quadNlim:", len(quadNlim[0]), len(quadNlim[1]), len(quadNlim[2]), len(quadNlim[3]), len(quadNlim[4]) 
	print "quadNrep:", len(quadNrep[0]), len(quadNrep[1]), len(quadNrep[2]), len(quadNrep[3]), len(quadNrep[4]) 
	print "quadnull:", len(quadnull[0]), len(quadnull[1]), len(quadnull[2]), len(quadnull[3]), len(quadnull[4]) 
	return allquad

def printSTDToFile(STD_Counts, NoutputFile, PoutputFile):
#Write STD	
#STD_all=[PupP_STD, PupN_STD, PdnP_STD, PdnN_STD, NupP_STD, NupN_STD, NdnP_STD, NdnN_STD]
	PupP=STD_Counts[0]
	PupN=STD_Counts[1]
	PdnP=STD_Counts[2]
	PdnN=STD_Counts[3]
	NupP=STD_Counts[4]
	NupN=STD_Counts[5]
	NdnP=STD_Counts[6]
	NdnN=STD_Counts[7]
	NSTD=dict(PupN.items() + PdnN.items() + NupN.items() + NdnN.items())
	PSTD=dict(PupP.items() + PdnP.items() + NupP.items() + NdnP.items())
	NWriter = open(NoutputFile, 'w')
	NWriter.write("Sequence")
	NWriter.write("\t")
	NWriter.write("S1")
	NWriter.write("\t")
	NWriter.write("S2")
	NWriter.write("\t")
	NWriter.write("S3")
	NWriter.write("\t")
	NWriter.write("S4")
	NWriter.write("\t")
	NWriter.write("S5")
	NWriter.write("\n")
	for key in NSTD: 
		NWriter.write(key)
		NWriter.write("\t")
		for item in NSTD[key]:
			NWriter.write(str(item))
			NWriter.write("\t")
		NWriter.write("\n")
	NWriter.close()
	
	PWriter = open(PoutputFile, 'w')
	PWriter.write("Sequence")
	PWriter.write("\t")
	PWriter.write("S1")
	PWriter.write("\t")
	PWriter.write("S2")
	PWriter.write("\t")
	PWriter.write("S3")
	PWriter.write("\t")
	PWriter.write("S4")
	PWriter.write("\t")
	PWriter.write("S5")
	PWriter.write("\n")
	for key in PSTD: 
		PWriter.write(key)
		PWriter.write("\t")
		for item in PSTD[key]:
			PWriter.write(str(item))
			PWriter.write("\t")
		PWriter.write("\n")
	PWriter.close()
	
def pickleIt(file_name,data, outdir):
	try:
		with open(outdir+file_name+".pickle", "wb") as output_file:
			pickle.dump(data, output_file,-1)
        	output_file.close()
	except Exception:
		print "Cannot open the file:",file_name
		
def unPickleAll(file_name, outdir):
	listHash=[]
	for name in file_name:
		dict={}
		dict=pickle.load(open(outdir+name+".pickle", "rb"))
		listHash.append(dict)
	return listHash 
	
def count_Quadrants2(STD_Counts, nN, nP): 
#	Input STD Counts list (broken down by Pup, Nup, etc, etc. 
#	STD_all=[PupP_STD, PupN_STD, PdnP_STD, PdnN_STD, NupP_STD, NupN_STD, NdnP_STD, NdnN_STD]
	#Compare by Pup/Pdn etc. etc. 
	Pup=STD_Counts[0:2]
	Pdn=STD_Counts[2:4]
	Nup=STD_Counts[4:6]
	Ndn=STD_Counts[6:8]
	All_Set=[Pup, Pdn, Nup, Ndn]
	
	quad1=[[],[],[],[],[],[]] # N and P limited (Co-limited) ; p>1 and n>1
	quad2=[[],[],[],[],[],[]] # P limited; p>1 and n<0
	quad3=[[],[],[],[],[], []] # N and P replete; p<0 and n<0
	quad4=[[],[],[],[],[],[]] # N limited 
	
	count=0
	for set in All_Set:
		for key in set[0]:
			count+=1
			Pnums=set[0][key]
			Nnums=set[1][key]
			for x in np.arange(5):
				p=Pnums[x]
				n=Nnums[x]
				if (p>=nP) & (n>=nN): 
					#Quadrant1 Co-limited
					quad1[x].append(key)
				elif (p>=nP) & (n<=nN):
					#Quadrant2 P-limited 
					quad2[x].append(key)
				elif (p<=nP) & (n<=nN):
					#Quadrant3 Replete
					quad3[x].append(key)
				elif (p<=nP) & (n>=nN):
					#Quadrant4 N-limited
					quad4[x].append(key)
	lenQuad1=[len(quad1[0]), len(quad1[1]), len(quad1[2]), len(quad1[3]), len(quad1[4])]
	lenQuad2=[len(quad2[0]), len(quad2[1]), len(quad2[2]), len(quad2[3]), len(quad2[4])] 
	lenQuad3=[len(quad3[0]), len(quad3[1]), len(quad3[2]), len(quad3[3]), len(quad3[4])]
	lenQuad4=[len(quad4[0]), len(quad4[1]), len(quad4[2]), len(quad4[3]), len(quad4[4])] 
	allLen=[lenQuad1, lenQuad2, lenQuad3, lenQuad4]
	allquad=[[quad1, quad2, quad3, quad4], allLen, count]
# 	print np.sum([len(quad1[0]),len(quad2[0]), len(quad3[0]), len(quad4[0]), len(quadPlim[0]), len(quadPrep[0]), len(quadNlim[0]), len(quadNrep[0]), len(quadnull[0])])							
# 	print "quad1:", len(quad1[0]), len(quad1[1]), len(quad1[2]), len(quad1[3]), len(quad1[4]) 
# 	print "quad2:", len(quad2[0]), len(quad2[1]), len(quad2[2]), len(quad2[3]), len(quad2[4]) 
# 	print "quad3:", len(quad3[0]), len(quad3[1]), len(quad3[2]), len(quad3[3]), len(quad3[4]) 
# 	print "quad4:", len(quad4[0]), len(quad4[1]), len(quad4[2]), len(quad4[3]), len(quad4[4]) 
	return allquad


def count_Quadrants3(STD_Counts, nN, nP): 
#	Input STD Counts list (broken down by Pup, Nup, etc, etc. 
#	STD_all=[PupP_STD, PupN_STD, PdnP_STD, PdnN_STD, NupP_STD, NupN_STD, NdnP_STD, NdnN_STD]
	#Compare by Pup/Pdn etc. etc. 
	Pup=STD_Counts[0:2]
	Pdn=STD_Counts[2:4]
	Nup=STD_Counts[4:6]
	Ndn=STD_Counts[6:8]
	All_Set=[Pup, Pdn, Nup, Ndn]
	
	quad1=[[],[],[],[],[],[]] # N and P limited (Co-limited) ; p>1 and n>1
	quad2=[[],[],[],[],[],[]] # P limited; p>1 and n<0
	quad3=[[],[],[],[],[],[]] # N and P replete; p<0 and n<0
	quad4=[[],[],[],[],[],[]] # N limited 
	
	count=0
	for set in All_Set:
		for key in set[0]:
			count+=1
			Pnums=set[0][key]
			Nnums=set[1][key]
			for x in np.arange(6):
				p=Pnums[x]
				n=Nnums[x]
				if (p>=nP) & (n>=nN): 
					#Quadrant1 Co-limited
					quad1[x].append(key)
				elif (p>=nP) & (n<=nN):
					#Quadrant2 P-limited 
					quad2[x].append(key)
				elif (p<=nP) & (n<=nN):
					#Quadrant3 Replete
					quad3[x].append(key)
				elif (p<=nP) & (n>=nN):
					#Quadrant4 N-limited
					quad4[x].append(key)
	lenQuad1=[len(quad1[0]), len(quad1[1]), len(quad1[2]), len(quad1[3]), len(quad1[4]), len(quad1[5])]
	lenQuad2=[len(quad2[0]), len(quad2[1]), len(quad2[2]), len(quad2[3]), len(quad2[4]), len(quad2[5])] 
	lenQuad3=[len(quad3[0]), len(quad3[1]), len(quad3[2]), len(quad3[3]), len(quad3[4]), len(quad3[5])]
	lenQuad4=[len(quad4[0]), len(quad4[1]), len(quad4[2]), len(quad4[3]), len(quad4[4]), len(quad4[5])] 
	allLen=[lenQuad1, lenQuad2, lenQuad3, lenQuad4]
	allquad=[[quad1, quad2, quad3, quad4], allLen, count]
# 	print np.sum([len(quad1[0]),len(quad2[0]), len(quad3[0]), len(quad4[0]), len(quadPlim[0]), len(quadPrep[0]), len(quadNlim[0]), len(quadNrep[0]), len(quadnull[0])])							
# 	print "quad1:", len(quad1[0]), len(quad1[1]), len(quad1[2]), len(quad1[3]), len(quad1[4]) 
# 	print "quad2:", len(quad2[0]), len(quad2[1]), len(quad2[2]), len(quad2[3]), len(quad2[4]) 
# 	print "quad3:", len(quad3[0]), len(quad3[1]), len(quad3[2]), len(quad3[3]), len(quad3[4]) 
# 	print "quad4:", len(quad4[0]), len(quad4[1]), len(quad4[2]), len(quad4[3]), len(quad4[4]) 
	return allquad


	