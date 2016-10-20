if(!exists('taxas'))source('analyzeNt.R')


humanFilter<-lapply(taxas,function(x)x[x[,'class']!='Mammalia'|is.na(x[,'class']),])
humanFilter<-humanFilter[!grepl('_Undetermined_',names(humanFilter))]
shortNames<-sub('_S[0-9]+_L[0-9]+_R.*$','',sub('blast_trim_','',names(humanFilter)))
taxaTab<-table(unlist(lapply(humanFilter,function(x)x$best)),rep(shortNames,sapply(humanFilter,nrow)))
propTab<-apply(taxaTab,2,function(x)x/sum(x))
cutTab<-propTab[apply(propTab,1,max)>.01,]
print(apply(taxaTab,2,sum))

breaks<-c(0,seq(min(cutTab)*.999,max(cutTab)*1.001,length.out=501))
cols<-c('white',tail(rev(heat.colors(520)),500))

pdf('out/heat.pdf',height=10,width=9)
heatmap(cutTab,col=cols,breaks=breaks,margins=c(5,9),scale='none')
dev.off()

insetLegend<-function(breaks,cols,insetPos=c(.015,.025,.25,.04)){
  insetPos<-c(grconvertX(insetPos[1],'nfc','user'),grconvertY(insetPos[2],'nfc','user'),grconvertX(insetPos[3],'nfc','user'),grconvertY(insetPos[4],'nfc','user'))
  breakPos<-((breaks[-1])-(min(breaks[-1])))/max((breaks[-1])-(min(breaks[-1])))*(insetPos[3]-insetPos[1])+insetPos[1]
  rect(breakPos[-1]+1e-3,insetPos[2],breakPos[-length(breakPos)],insetPos[4],col=cols[-1],xpd=NA,border=NA)
  rect(insetPos[1],insetPos[2],insetPos[3],insetPos[4],xpd=NA)
  prettyLabs<-pretty(breaks)
  prettyLabs<-prettyLabs[prettyLabs<max(breaks)]
  prettyPos<-prettyLabs
  prettyPos<-(prettyLabs-(min(breaks[-1])))/((max(breaks[-1]))-(min(breaks[-1])))*(insetPos[3]-insetPos[1])+insetPos[1]
  segments(prettyPos,insetPos[2],prettyPos,insetPos[2]-diff(insetPos[c(2,4)])*.1,xpd=NA)
  text(prettyPos,insetPos[2]-diff(insetPos[c(2,4)])*.175,prettyLabs,xpd=NA,adj=c(.5,1),cex=.85)
  text(mean(insetPos[c(1,3)]),insetPos[4]+diff(insetPos[c(2,4)])*.45,"Proportion of sample",xpd=NA,adj=c(.5,0))
}

types<-sapply(strsplit(colnames(cutTab),'-'),'[',1)
types2<-sapply(strsplit(colnames(cutTab),'-'),'[',3)
newOrder<-order(types,ifelse(is.na(types2),'ZZ',types2))
cutTab<-cutTab[order(apply(cutTab,1,sum)),newOrder]
types<-types[newOrder]
types2<-types2[newOrder]
notInControl<-which(apply(cutTab[,types=='h20'],1,sum)==0 &apply(cutTab[,types!='h20']>0,1,sum)>1)
#cutTab<-cutTab[order(apply(cutTab,1,sum)),]
pdf('out/heatSort.pdf',width=9,height=16)
  par(mar=c(.1,15,7,.1))
  image(1:ncol(cutTab),1:nrow(cutTab),t(cutTab),col=cols,breaks=breaks,xlab='',ylab='',xaxt='n',yaxt='n')
  axis(3,1:ncol(cutTab),colnames(cutTab),las=2)
  axis(2,1:nrow(cutTab),rownames(cutTab),las=1)
  abline(v=which(types[-1]!=types[-length(types)])+.5)
  abline(h=2:nrow(cutTab)-.5,col='#00000033')
  abline(h=2:nrow(cutTab)-.5,col='#00000033')
  abline(v=2:ncol(cutTab)-.5,col='#00000011')
  points(rep(convertLineToUser(.75,2),length(notInControl)),notInControl,pch='*',xpd=NA,cex=2,col='red')
  box()
  insetLegend(breaks,cols,insetPos=c(.05,.93,.25,.945))
dev.off()

virusHits<-mapply(function(tax,reads)tax[tax$qName %in% reads,],taxonomy,lapply(taxas,function(x)rownames(x)[x$best=='Viruses']),SIMPLIFY=FALSE)

