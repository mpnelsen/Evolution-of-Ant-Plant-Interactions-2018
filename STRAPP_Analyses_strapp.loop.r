require(BAMMtools)
require(ape)

mcmcout <- read.csv(file="/PATH/mcmc_out_medusa.txt", header=T)
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

master<-read.csv(file="/PATHTODATAFILE",stringsAsFactors=FALSE)
tree<-read.tree("/PATHTOTREE")
edata <- getEventData(tree, eventdata = "/PATH/event_data_medusa.txt", burnin=0.1,nsamples=2000)


master<-master[master$Family %in% "Formicidae",]
master<-master[,c("Taxon","DietModanyplant","ForagingModallupper","NestingModallupper")]
master[master=="?"]<-NA
master[,2]<-as.numeric(master[,2])
master[,3]<-as.numeric(master[,3])
master[,4]<-as.numeric(master[,4])
master$CombinedTrait<-NA
master$DietArbOnly<-NA
master$CombinedTrait[master$DietModanyplant %in% 0 & master$ForagingModallupper %in% 0 & master$NestingModallupper %in% 0]<-1
master$CombinedTrait[master$DietModanyplant %in% 1 & master$ForagingModallupper %in% 0 & master$NestingModallupper %in% 0]<-2
master$CombinedTrait[master$DietModanyplant %in% 0 & master$ForagingModallupper %in% 1 & master$NestingModallupper %in% 0]<-3
master$CombinedTrait[master$DietModanyplant %in% 0 & master$ForagingModallupper %in% 0 & master$NestingModallupper %in% 1]<-4
master$CombinedTrait[master$DietModanyplant %in% 1 & master$ForagingModallupper %in% 1 & master$NestingModallupper %in% 0]<-5
master$CombinedTrait[master$DietModanyplant %in% 1 & master$ForagingModallupper %in% 0 & master$NestingModallupper %in% 1]<-6
master$CombinedTrait[master$DietModanyplant %in% 0 & master$ForagingModallupper %in% 1 & master$NestingModallupper %in% 1]<-7
master$CombinedTrait[master$DietModanyplant %in% 1 & master$ForagingModallupper %in% 1 & master$NestingModallupper %in% 1]<-8
master$CombinedTrait[rowSums(is.na(master[,c("DietModanyplant","ForagingModallupper","NestingModallupper")]))>0]<-NA

master$DietArbOnly[master$CombinedTrait %in% 1]<-1
master$DietArbOnly[master$CombinedTrait %in% 2]<-2
master$DietArbOnly[master$CombinedTrait %in% c(3,4,7)]<-3
master$DietArbOnly[master$CombinedTrait %in% c(5,6,8)]<-4
master$DietArbOnly[is.na(master$CombinedTrait)]<-NA

mods<-c("CombinedTrait","DietArbOnly")

df<-as.data.frame(matrix(nrow=(length(mods)*3),ncol=13))
colnames(df)<-c("CodingScheme","Estimate1","Estimate2","Estimate3","Estimate4","Estimate5","Estimate6","Estimate7","Estimate8","Pval","Method","TwoTailed","Rate")

for(x in 1:length(mods)){
	tr.vec<-master[,mods[x]]
	names(tr.vec)<-master$Taxon
	keeps<-names(tr.vec)[tr.vec %in% c(1,2,3,4,5,6,7,8)]
	#tr.vec[tr.vec %in% "?"]<-NA
	tr.vec<-tr.vec[names(tr.vec) %in% keeps]
	edata.sub<-subtreeBAMM(edata,tips=keeps)
	rates<-c("net diversification","speciation","extinction")
	#when netdiv, may need to use logrates=FALSE
	if(any(!unique(tr.vec) %in% c(0,1))){
		if(length(unique(tr.vec))==4){
			#when 3 states or more use method="kruskal"
			#need to make factor
			tr.vec<-factor(tr.vec)
			for(i in 1:length(rates)){
				strapp.res<-traitDependentBAMM(ephy=edata.sub,traits=tr.vec,reps=1000,logrates=FALSE,method="kruskal",rate=rates[i])
				df[(x*3)-3+i,1]<-paste(mods[x],rates[i],sep="_")
				df[(x*3)-3+i,1]<-gsub(" ","",df[(x*3)-3+i,1])
				df[(x*3)-3+i,2]<-strapp.res$estimate$'1'
				df[(x*3)-3+i,3]<-strapp.res$estimate$'2'
				df[(x*3)-3+i,4]<-strapp.res$estimate$'3'
				df[(x*3)-3+i,5]<-strapp.res$estimate$'4'
				df[(x*3)-3+i,10]<-strapp.res$p.val
				df[(x*3)-3+i,11]<-strapp.res$method
				df[(x*3)-3+i,12]<-strapp.res$two.tailed
				df[(x*3)-3+i,13]<-strapp.res$rate			
			}
		}
		if(length(unique(tr.vec))==8){
			#when 3 states or more use method="kruskal"
			#need to make factor
			tr.vec<-factor(tr.vec)
			for(i in 1:length(rates)){
				strapp.res<-traitDependentBAMM(ephy=edata.sub,traits=tr.vec,reps=1000,logrates=FALSE,method="kruskal",rate=rates[i])
				df[(x*3)-3+i,1]<-paste(mods[x],rates[i],sep="_")
				df[(x*3)-3+i,1]<-gsub(" ","",df[(x*3)-3+i,1])
				df[(x*3)-3+i,2]<-strapp.res$estimate$'1'
				df[(x*3)-3+i,3]<-strapp.res$estimate$'2'
				df[(x*3)-3+i,4]<-strapp.res$estimate$'3'
				df[(x*3)-3+i,5]<-strapp.res$estimate$'4'
				df[(x*3)-3+i,6]<-strapp.res$estimate$'5'
				df[(x*3)-3+i,7]<-strapp.res$estimate$'6'
				df[(x*3)-3+i,8]<-strapp.res$estimate$'7'
				df[(x*3)-3+i,9]<-strapp.res$estimate$'8'
				df[(x*3)-3+i,10]<-strapp.res$p.val
				df[(x*3)-3+i,11]<-strapp.res$method
				df[(x*3)-3+i,12]<-strapp.res$two.tailed
				df[(x*3)-3+i,13]<-strapp.res$rate
			}
		}
	}
	if(all(unique(tr.vec) %in% c(0,1))){
		#when 2 states or more use method="mann-whitney"
		for(i in 1:length(rates)){
			strapp.res<-traitDependentBAMM(ephy=edata.sub,traits=tr.vec,reps=1000,logrates=FALSE,method="mann-whitney",rate=rates[i])
			df[(x*3)-3+i,1]<-paste(mods[x],rates[i],sep="_")
			df[(x*3)-3+i,1]<-gsub(" ","",df[(x*3)-3+i,1])
			df[(x*3)-3+i,2]<-strapp.res$estimate$'0'
			df[(x*3)-3+i,3]<-strapp.res$estimate$'1'
			df[(x*3)-3+i,5]<-strapp.res$p.val
			df[(x*3)-3+i,6]<-strapp.res$method
			df[(x*3)-3+i,7]<-strapp.res$two.tailed
			df[(x*3)-3+i,8]<-strapp.res$rate
		}
	}
}
any(df$Pval<0.05)
write.csv(df,file="/PATH/summmary.strapp.analyses.csv",row.names=FALSE)
