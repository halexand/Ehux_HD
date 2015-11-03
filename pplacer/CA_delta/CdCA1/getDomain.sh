##Quick bash foo script to get a subset of fasta from hmm search etc. 
##Input requires a tab delimited file QueryID, Start, Stop
## bash getDomain.sh < inputFile
database=~/Analysis/MMETSP_BlastDB/MMETSP_150604_aminoacid/MMETSP_150604
database=../../../Blastdb/EmiHux_Clean_Proteins

while read entry start stop
    do 
        start_adj=$(($start+1))
        blastdbcmd -db $database -entry $entry -range $start_adj-$stop
done 
