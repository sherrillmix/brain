
library(dnar)
library(taxonomizr)
library(parallel)

getNamesAndNodes()
taxaNodes<-read.nodes('nodes.dmp')
taxaNames<-read.names('names.dmp')
if(!file.exists('accessionTaxa.sql')){
  tmp<-tempdir()
  sqlFile<-'accessionTaxa.sql'
  #this is a big download
  getAccession2taxid(tmp)
  read.accession2taxid(list.files(tmp,'accession2taxid.gz$',full.names=TRUE),sqlFile)
  file.remove(list.files(tmp,'accession2taxid.gz$',full.names=TRUE))
}

blastFiles<-list.files('work','.blast.gz$',full.names=TRUE)
taxas<-mclapply(blastFiles,function(ii){
  message(ii)
  outFile<-sub('.blast.gz$','_taxa.csv',ii)
  if(file.exists(outFile)){
    message(' Reading ',outFile)
    taxaAssigns<-read.csv(outFile,row.names=1,stringsAsFactors=FALSE)
  }else{
    message(' Creating ',outFile)
    x<-read.blast(ii)
    x$accession<-sapply(strsplit(x$tName,'\\|'),'[[',4)
    x$taxa<-accessionToTaxa(x$accession,sqlFile)
    x$maxBit<-ave(x$score,x$qName,FUN=max)
    x<-x[x$score==x$maxBit&!is.na(x$taxa),]
    taxonomy<-getTaxonomy(x$taxa,taxaNodes,taxaNames,mc.cores=1)
    taxaAssigns<-do.call(rbind,by(as.data.frame(taxonomy,stringsAsFactors=FALSE),x$qName,FUN=condenseTaxa))
    taxaAssigns$best<-apply(taxaAssigns,1,lastNotNa)
    bestScore<-x[!duplicated(x$qName),c('qName','alignLength','percID')]
    rownames(bestScore)<-bestScore$qName
    taxaAssigns<-cbind(taxaAssigns,bestScore[rownames(taxaAssigns),c('alignLength','percIdent')])
    write.csv(taxaAssigns,outFile)
  }
  return(taxaAssigns)
},mc.cores=15)
names(taxas)<-basename(blastFiles)

