#Adds outgroup info to spreadsheet and also updates tribe taxonomy - if tribe is listed as incertae sedis, the genus is assigned to be the tribe.

master<-read.csv(file="/PATH/combined.updated.updated.pruned.csv",stringsAsFactors=FALSE)
outs<-read.csv(file="/PATH/standards.w.taxonomy.dupaccsrem.alltax.wouts.outsonly.csv",stringsAsFactors=FALSE)
outs<-outs[,c(1:10,12:18)]
combo<-merge(master,outs,all=TRUE)
for(v in 1:nrow(combo)){
	if(combo$Tribe[v]=="Incertae_Sedis"){
		combo$Tribe[v]<-combo$Genus[v]
	}
}
combo$Order<-"Hymenoptera"
write.csv(combo,file="/PATH/combined.updated.updated.pruned.wouts.csv",row.names=FALSE)
