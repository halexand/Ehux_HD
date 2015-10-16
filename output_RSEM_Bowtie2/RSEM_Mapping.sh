wd=/N/dc2/projects/ldeo/Ehux_HD/rsem_Bowtie2

cd $wd

forward=E1-Q_ATCACG__trimmed_paired_1.fastq
reverse=${forward/1.fastq/2.fastq}

n=${forward%__trimmed_paired_1.fastq}

RSEM-1.2.20/rsem-calculate-expression --paired-end \
   -p 8 \
   --bowtie2 \
   --bowtie2-mismatch-rate 0.2 $forward \
   $reverse \
   Ehux2 \
   ${n}


