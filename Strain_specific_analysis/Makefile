all: Species_GeneLists_Venn 

clean: 
	rm *	

Species_GeneLists_Venn: Ehux_All_Transcripts_Cleaned_RSEM.nt.fa.RSEM.tab GetVenn_Species.py
	mkdir Species_GeneLists_Venn
	python GetVenn_Species.py
