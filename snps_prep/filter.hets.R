#!/usr/bin/Rscript

data<-read.table("out.hwe",header=T)
split <- data.frame(do.call('rbind', strsplit(as.character(data$OBS.HOM1.HET.HOM2.),'/',fixed=TRUE)))
part <-data[,1:2]
ready <-cbind(part,split)
names(ready)<-c("chr","pos","homo1","het","homo2")
ready$homo1<-as.numeric(as.character(ready$homo1))
ready$homo2<-as.numeric(as.character(ready$homo2))
ready$het<-as.numeric(as.character(ready$het))
ready$obs_het<-0
for (i in nrow(ready)){
  ready$obs_het<-ready$het/(ready$homo1+ready$het+ready$homo2)
}
hist(ready$obs_het,breaks=60)
na<-subset(ready,ready$obs_het=="NaN")
het<-subset(ready,ready$obs_het>=0.6)
remove<-rbind(na,het)
remove_2col<-remove[,1:2]
write.table(remove_2col,file="snps_het06_filter",sep="\t",quote=FALSE,row.names=F)