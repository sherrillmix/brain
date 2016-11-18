
library(dnar)
library(taxonomizr)
library(parallel)

getNamesAndNodes()
taxaNodes<-read.nodes('nodes.dmp')
taxaNames<-read.names('names.dmp')
sqlFile<-'accessionTaxa.sql'
if(!file.exists(sqlFile)){
  tmp<-tempdir()
  #this is a big download
  getAccession2taxid(tmp)
  read.accession2taxid(list.files(tmp,'accession2taxid.gz$',full.names=TRUE),sqlFile)
  file.remove(list.files(tmp,'accession2taxid.gz$',full.names=TRUE))
}

blastFiles<-list.files('work','.blast.gz$',full.names=TRUE)
taxas<-lapply(blastFiles,function(ii){
  message(ii)
  outFile<-sub('.blast.gz$','_taxa.csv',ii)
  outFile2<-sub('.blast.gz$','_allHits.csv',ii)
  if(file.exists(outFile)&file.exists(outFile2)){
    message(' Reading ',outFile)
    taxaAssigns<-read.csv(outFile,row.names=1,stringsAsFactors=FALSE)
    taxonomy<-read.csv(outFile2,stringsAsFactors=FALSE)
  }else{
    message(' Creating ',outFile)
    message('  Reading blast')
    x<-read.blast(ii)
    x$accession<-sapply(strsplit(x$tName,'\\|'),'[[',4)
    message('  Accession to taxonomy')
    x$taxa<-accessionToTaxa(x$accession,sqlFile)
    x$maxBit<-ave(x$score,x$qName,FUN=max)
    x<-x[x$score==x$maxBit&!is.na(x$taxa),]
    gc()
    message('  Getting upstream taxonomy')
    taxonomy<-getTaxonomy(x$taxa,taxaNodes,taxaNames,mc.cores=10)
    taxonomy<-as.data.frame(taxonomy,stringsAsFactors=FALSE)
    message('  Condensing taxonomy')
    taxaAssigns<-do.call(rbind,by(taxonomy,x$qName,FUN=condenseTaxa))
    taxonomy$qName<-x$qName
    write.csv(taxonomy,outFile2)
    taxaAssigns$best<-apply(taxaAssigns,1,lastNotNa)
    bestScore<-x[!duplicated(x$qName),c('qName','alignLength','percID')]
    rownames(bestScore)<-bestScore$qName
    taxaAssigns<-cbind(taxaAssigns,bestScore[rownames(taxaAssigns),c('alignLength','percID')])
    write.csv(taxaAssigns,outFile)
  }
  return(list('taxa'=taxaAssigns,'taxonomy'=taxonomy))
})

taxonomy<-lapply(taxas,'[[','taxonomy')
taxas<-lapply(taxas,'[[','taxa')
names(taxas)<-names(taxonomy)<-basename(blastFiles)
