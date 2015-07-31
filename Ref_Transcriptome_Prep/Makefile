all: cosci.x.diatoms.xml diatoms.x.cosci.xml cosci.x.diatoms.besthit.homol diatoms.x.cosci.besthit.orthol cosci.x.diatom.reciprocal.besthits.ortho
# Make file for the initial annotation of Cosci transcriptome based on other diatom genomes from the adapted khmer protocol
# 
#
clean: 
	rm *.xml 

cosci.x.diatoms.xml: 
	cp ../cosci.x.diatoms.xml . 

diatoms.x.cosci.xml:
	cp ../diatoms.x.cosci.xml .

cosci.x.diatoms.besthit.homol: cosci.x.diatoms.xml Blast_besthit.cosci.py
	python Blast_besthit.cosci.py cosci.x.diatoms.xml cosci.x.diatoms.besthit.homol	
#Get the best blast hits for cosci against the diatoms; save as a pickled hash
diatoms.x.cosci.besthit.orthol: diatoms.x.cosci.xml Blast_besthit.diatom.py
	python Blast_besthit.diatom.py diatoms.x.cosci.xml diatoms.x.cosci.besthit.orthol
#Get the best blast hits for diatoms against cosci; save as a pickled hash

cosci.x.diatom.reciprocal.besthits.ortho: Parse_reciprocal_besthits.py cosci.x.diatoms.besthit.homol diatoms.x.cosci.besthit.orthol
	python Parse_reciprocal_besthits.py cosci.x.diatoms.besthit.homol diatoms.x.cosci.besthit.orthol cosci.x.diatom.reciprocal.besthits.ortho
#Create a list of orthologs for the cosci annotation

cosci.x.diatoms.besthit.homol.tab: BestHit_Homologues_toTab.py cosci.x.diatoms.besthit.homol diatoms.x.cosci.besthit.orthol
	python BestHit_Homologues_toTab.py
