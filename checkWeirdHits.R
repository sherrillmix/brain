if(!exists('taxas'))source('analyzeNt.R')
library(dnaplotr)
library(bio3d)

virusHits<-mapply(function(tax,reads)tax[tax$qName %in% reads,],taxonomy,lapply(taxas,function(x)rownames(x)[x$best=='Viruses']),SIMPLIFY=FALSE)
ralstoniaHits<-mapply(function(tax,reads)tax[tax$qName %in% reads,],taxonomy,lapply(taxas,function(x)rownames(x)[x$best=='Ralstonia solanacearum']),SIMPLIFY=FALSE)

trims<-sub('blast.gz$','fastq.gz',sub('^blast_','',names(taxas)))
bases<-sub('_R[12][_.].*','',trims)

ralReads<-mclapply(unique(bases),function(thisBase){
  message(thisBase)
  baseSelect<-bases==thisBase
  fastqs<-file.path('work',trims[baseSelect])
  targets<-unique(unlist(sapply(ralstoniaHits[baseSelect],function(xx)xx$qName)))
  reads<-do.call(cbind,lapply(fastqs,function(xx,targets){
    out<-read.fastq(xx)
    rownames(out)<-out$name
    out[targets,c('seq','qual')]
  },targets))
  return(reads)
},mc.cores=2)
names(ralReads)<-unique(bases)

withAs(xx=ralReads[['trim_Phi6Spikec']],write.fastq(rownames(xx),xx[,1],xx[,2],'out/phi6_ral1.fastq.gz'))
withAs(xx=ralReads[['trim_Phi6Spikec']],write.fastq(rownames(xx),xx[,3],xx[,4],'out/phi6_ral1.fastq.gz'))

aligned<-apply(seqaln(seqSplit(ralReads[['trim_Phi6Spikec']][sample(nrow(ralReads[['trim_Phi6Spikec']]),2000),1],fill='-'),protein=FALSE,extra.args='')$ali,1,paste,collapse='')
aligned2<-apply(seqaln(seqSplit(ralReads[['trim_test']][,1],fill='-'),protein=FALSE,extra.args='')$ali,1,paste,collapse='')
png('out/phi6_ralstonia.png',width=2000,height=2000,res=250);plotDNA(sort(aligned));dev.off()
png('out/test_ralstonia.png',width=2000,height=2000,res=250);plotDNA(sort(aligned2));dev.off()
