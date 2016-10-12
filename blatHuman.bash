#!/bin/bash
set -e
date
for ii in data/*R[12]*.fastq.gz;do
  echo $ii 
  base=`basename $ii`
  out=work/${base%.fastq.gz}.blat.gz
  if [ ! -e "$out" ];then
    echo Running blat on $ii $(date)
    sem --no-notice -j15 "zcat $ii|awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,\">\");print}; if(P==4)P=0; P++}'|~/installs/blat/blat ~/installs/blat/hg38.2bit stdin >(gzip>$out)"
    echo Done with blat on $ii $(date)
  else
    echo File already exists. Skipping
  fi
done
sem --wait
date
