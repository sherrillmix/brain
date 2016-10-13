library(parallel)
library(dnar)
blatFiles<-list.files('work','blat.gz$',full.names=TRUE)
filters<-mclapply(blatFiles,function(file){
  blat<-read.blat(file)
  return(unique(blat$qName[blat$match/blat$qSize>.5]))
},mc.cores=15)
names(filters)<-blatFiles
baseNames<-sub('_R[12]_.*$','',blatFiles)

mclapply(unique(baseNames),function(base){
  #from dnapy package
  files<-blatFiles[baseNames==base]
  filter<-unique(unlist(filters[files]))
  tmp<-tempfile()
  writeLines(filter,tmp)
  datFiles<-sub('blat.gz','fastq.gz',sub('work/','data/',files))
  outFiles<-sub('data/','work/trim_',datFiles)
  cmd<-sprintf('removereads -f %s %s -o %s',tmp,paste(datFiles,collapse=' '),paste(outFiles,collapse=','))
  message(cmd)
  system(cmd)
  file.remove(tmp)
},mc.cores=10)
