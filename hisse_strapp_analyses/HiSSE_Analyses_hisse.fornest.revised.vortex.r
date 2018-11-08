require(ape)
require(phytools)
require(hisse)
require(parallel)
#slightly modified from SI code in Beaulieu & O'Meara (2016)

	
SimFunction <- function(phy=NULL,sim.dat=NULL,sampling.file=NULL,file.name=NULL){
	sampling.file[,2]<-sampling.file[,1]
	sampling.file.new <- sampling.file[phy$tip.label,]
	sampling.matrix <- as.matrix(sampling.file.new[,1])
	sampling.f <- as.vector(sampling.matrix)
	trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
	trans.rates.hisse <- ParDrop(trans.rates.hisse, c(3,5,8,10))
	trans.rates.hisse.test <- TransMatMaker(hidden.states=TRUE)
	trans.rates.hisse.test <- ParDrop(trans.rates.hisse.test, c(3,5,8,9,10,12))
	trans.rates.hisse.red <- trans.rates.hisse
	trans.rates.hisse.red.test <- trans.rates.hisse.test
	trans.rates.hisse.red[!is.na(trans.rates.hisse.red) & !trans.rates.hisse.red == 0] = 1
	trans.rates.hisse.red.test[!is.na(trans.rates.hisse.red.test) & !trans.rates.hisse.red.test == 0] = 1
#	trans.rates.hisse.red.test[2,1]=2
#	trans.rates.hisse.red.test[1,2]=3
	trans.rates.bisse <- TransMatMaker(hidden.states=FALSE)
	trans.rates.bisse.red <- trans.rates.bisse
	trans.rates.bisse.red[!is.na(trans.rates.bisse.red)] = 1
	
	hisse.fit <- NA
	bisse.fit <- NA
	
	RunModel <- function(model.number, sampling.f){
		if(model.number==1){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.bisse))	
		}
		if(model.number==2){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.bisse))
		}
		if(model.number==3){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.bisse.red))	
		}
		if(model.number==4){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.bisse.red))
		}
		if(model.number==5){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,2,2), trans.rate=trans.rates.hisse.red))
		}
		if(model.number==6){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red))
		}
		if(model.number==7){
			try(hisse.fit <- hisse.null4(phy, sim.dat, f=sampling.f, bounded.search=TRUE))
		}
		if(model.number==8){
			try(hisse.fit <- hisse.null4(phy, sim.dat, f=sampling.f, bounded.search=TRUE, eps.anc=rep(1,8)))
		}
		if(model.number==9){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.hisse.red))
		}
		if(model.number==10){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red))
		}
		if(model.number==11){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,2), trans.rate=trans.rates.hisse.red))
		}
		if(model.number==12){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red))		
		}
		if(model.number==13){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.hisse.red.test))
		}
		if(model.number==14){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red.test))
		}
		if(model.number==15){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,2), trans.rate=trans.rates.hisse.red.test))
		}
		if(model.number==16){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red.test))	
		}
		if(model.number==17){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,2,1,3), eps.anc=c(1,2,1,3), trans.rate=trans.rates.hisse.red))
		}
		if(model.number==18){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,2,1,3), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red))	
		}
		if(model.number==19){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,2,1,3), eps.anc=c(1,2,1,3), trans.rate=trans.rates.hisse.red.test))
		}
		if(model.number==20){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,2,1,3), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red.test))	
		}
		if(model.number==21){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,2,3), trans.rate=trans.rates.hisse.red))
		}
		if(model.number==22){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red))	
		}
		if(model.number==23){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,2,3), trans.rate=trans.rates.hisse.red.test))
		}
		if(model.number==24){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red.test))	
		}		
		save(phy, hisse.fit, file=paste(file.name, model.number, "Rsave", sep="."))
	}
	mclapply(1:24, RunModel, sampling.f=sampling.f, mc.cores=12)

}


tree<-read.tree("TREEFILE")

spec<-read.csv(file="TRAITDATAFILE",stringsAsFactors=FALSE)

sampling.fractions<-read.delim("SAMPLINGFRACTIONSFILE")
sampling.fractions<-data.frame(sampling.fractions[,2],row.names=sampling.fractions[,1],stringsAsFactors=FALSE)
colnames(sampling.fractions)<-"f"
datasets<-c("DietModslashboth","DietModslashlower","DietModslashupper","DietModalllower","DietModallupper","DietModanyplant","DietModonlyplant","DietModonlypred", "DietModanypred", "DietModonlyomni", "DietModanyomni","ForagingMod","ForagingModalllower","ForagingModallupper","ForagingModanycanopy","ForagingModanyground","NestingMod","NestingModalllower","NestingModallupper","NestingModanycanopy","NestingModanyground")
diet.datasets<-c("DietModslashboth","DietModslashlower","DietModslashupper","DietModalllower","DietModallupper","DietModanyplant","DietModonlyplant","DietModonlypred", "DietModanypred", "DietModonlyomni", "DietModanyomni")
diet.datasets.nopoly.first<-c("DietModalllower","DietModallupper","DietModanyplant","DietModonlyplant")
diet.datasets.nopoly.second<-c("DietModonlypred", "DietModanypred", "DietModonlyomni", "DietModanyomni")
fornest.datasets<-c("ForagingMod","ForagingModalllower","ForagingModallupper","ForagingModanycanopy","NestingMod","NestingModalllower","NestingModallupper","NestingModanycanopy","ForagingModalllower","NestingModalllower","ForagingModanyground","NestingModanyground")
fornest.datasets.nopoly<-c("ForagingModanycanopy","ForagingModanyground","NestingModanycanopy","NestingModanyground")

for(pq in 1:length(fornest.datasets.nopoly)){
	df<-data.frame(spec[,fornest.datasets.nopoly[pq]],row.names=spec$Taxon,stringsAsFactors=FALSE)
	colnames(df)<-fornest.datasets.nopoly[pq]
	drops<-rownames(df)[df[,1]=="?"]
	df.wo.drops<-data.frame(rownames(df)[df[,1]!="?"],df[df[,1]!="?",],row.names=rownames(df)[df[,1]!="?"],stringsAsFactors=FALSE)
	colnames(df.wo.drops)<-c("Genus_species","states")
	tree.wo.drops<-drop.tip(tree,drops)
	sampling.fractions.mod<-data.frame(sampling.fractions[!rownames(sampling.fractions) %in% drops,],row.names=rownames(sampling.fractions)[!rownames(sampling.fractions) %in% drops],stringsAsFactors=FALSE)
	colnames(sampling.fractions.mod)<-"f"
	SimFunction(phy=tree.wo.drops,sim.dat=df.wo.drops,sampling.file=sampling.fractions.mod,file.name=paste("/PATH/",fornest.datasets.nopoly[pq],".ex.model",sep=""))
	EmpiricalRunSummary <- function(file.name){
		models <- system(paste("ls -1 ", file.name, "*.Rsave", sep=""), intern=TRUE)
		for(i in 1:length(models)){
			load(models[i])
			write.table(t(c(hisse.fit$loglik, hisse.fit$AIC, max(hisse.fit$index.par)-1, hisse.fit$solution)), file=paste("/PATH/",fornest.datasets.nopoly[pq],"_empirical.results",sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE)
		}
	}
	EmpiricalRunSummary(paste("/PATH/",fornest.datasets.nopoly[pq],".ex.model",sep=""))
}

for(pq in 1:length(diet.datasets.nopoly.first)){
	df<-data.frame(spec[,diet.datasets.nopoly.first[pq]],row.names=spec$Taxon,stringsAsFactors=FALSE)
	colnames(df)<-diet.datasets.nopoly.first[pq]
	drops<-rownames(df)[df[,1]=="?"]
	df.wo.drops<-data.frame(rownames(df)[df[,1]!="?"],df[df[,1]!="?",],row.names=rownames(df)[df[,1]!="?"],stringsAsFactors=FALSE)
	colnames(df.wo.drops)<-c("Genus_species","states")
	tree.wo.drops<-drop.tip(tree,drops)
	sampling.fractions.mod<-data.frame(sampling.fractions[!rownames(sampling.fractions) %in% drops,],row.names=rownames(sampling.fractions)[!rownames(sampling.fractions) %in% drops],stringsAsFactors=FALSE)
	colnames(sampling.fractions.mod)<-"f"
	SimFunction(phy=tree.wo.drops,sim.dat=df.wo.drops,sampling.file=sampling.fractions.mod,file.name=paste("/PATH/",diet.datasets.nopoly.first[pq],".ex.model",sep=""))
	EmpiricalRunSummary <- function(file.name){
		models <- system(paste("ls -1 ", file.name, "*.Rsave", sep=""), intern=TRUE)
		for(i in 1:length(models)){
			load(models[i])
			write.table(t(c(hisse.fit$loglik, hisse.fit$AIC, max(hisse.fit$index.par)-1, hisse.fit$solution)), file=paste("/PATH/",diet.datasets.nopoly.first[pq],"_empirical.results",sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE)
		}
	}
	EmpiricalRunSummary(paste("/PATH/",diet.datasets.nopoly.first[pq],".ex.model",sep=""))
}

for(pq in 1:length(diet.datasets.nopoly.second)){
	df<-data.frame(spec[,diet.datasets.nopoly.second[pq]],row.names=spec$Taxon,stringsAsFactors=FALSE)
	colnames(df)<-diet.datasets.nopoly.second[pq]
	drops<-rownames(df)[df[,1]=="?"]
	df.wo.drops<-data.frame(rownames(df)[df[,1]!="?"],df[df[,1]!="?",],row.names=rownames(df)[df[,1]!="?"],stringsAsFactors=FALSE)
	colnames(df.wo.drops)<-c("Genus_species","states")
	tree.wo.drops<-drop.tip(tree,drops)
	sampling.fractions.mod<-data.frame(sampling.fractions[!rownames(sampling.fractions) %in% drops,],row.names=rownames(sampling.fractions)[!rownames(sampling.fractions) %in% drops],stringsAsFactors=FALSE)
	colnames(sampling.fractions.mod)<-"f"
	SimFunction(phy=tree.wo.drops,sim.dat=df.wo.drops,sampling.file=sampling.fractions.mod,file.name=paste("/PATH/",diet.datasets.nopoly.second[pq],".ex.model",sep=""))
	EmpiricalRunSummary <- function(file.name){
		models <- system(paste("ls -1 ", file.name, "*.Rsave", sep=""), intern=TRUE)
		for(i in 1:length(models)){
			load(models[i])
			write.table(t(c(hisse.fit$loglik, hisse.fit$AIC, max(hisse.fit$index.par)-1, hisse.fit$solution)), file=paste("/PATH/",diet.datasets.nopoly.second[pq],"_empirical.results",sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE)
		}
	}
	EmpiricalRunSummary(paste("/PATH/",diet.datasets.nopoly.second[pq],".ex.model",sep=""))
}
