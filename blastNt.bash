#!/bin/bash
set -e
date
for ii in work/trim*R[12]*.fastq.gz;do
  echo $ii 
  base=`basename $ii`
  out=work/blast_${base%.fastq.gz}.blast.gz
  if [ ! -e "$out" ];then
    echo Running blast on $ii $(date)
    zcat $ii|awk 'BEGIN{P=1}{if(P==1||P==2){gsub("^@",">");print}; if(P==4)P=0; P++}'|parallel --block 100k  -L2 -j4 --pipe blastn -db ~/db/nt/nt -outfmt 6 -num_threads 5 -query - |gzip>$out
  else
    echo File already exists. Skipping
  fi
done
date
