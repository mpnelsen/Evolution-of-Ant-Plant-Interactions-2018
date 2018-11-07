#This then edits the retrieved file and compares taxon names obtained from NCBI with those in a list derived from AntCat.  
#Taxonomic discrepancies between the NCBI file and AntCat file were checked and names updated.
#It then creates a separate file for each locus containing the accession numbers to retrieve from NCBI.

#This checks which NCBI taxa retrieved are not in the AntCat list of taxa with geographic data
new<-read.csv(file="Formicidae.combined.acc.nos.csv",stringsAsFactors=FALSE)
new$Taxon<-gsub(" ","_",new$Taxon)
old<-read.csv(file="antcat_3june2016_clean2_geography.csv",stringsAsFactors=FALSE,row.names=1)
new.not.old<-new$Taxon[!new$Taxon %in% old$Taxon]
length(new.not.old)
write.csv(new.not.old,file="Formicidae.combined.acc.nos.ABSENT.OLD.csv")
##########
##########
yesterday<-read.csv(file="Formicidae.combined.acc.nos.ABSENT.OLD.updated.csv",stringsAsFactors=FALSE,row.names=1,na.strings=c("","NA"))

#this file was then manually modified to record which taxa not in the antcat geo checklist to keep or drop
#it also includes conversions between old and new names and links the two lists
today<-read.csv(file="Formicidae.combined.acc.nos.ABSENT.OLD.csv",stringsAsFactors=FALSE,row.names=1,na.strings=c("","NA"))
colnames(today)<-"Old"
nrow(yesterday)
nrow(today)
sum(!today$Old %in% yesterday$Old)
today$Old[!today$Old %in% yesterday$Old]
today$Current<-NA
today$Notes<-NA
today$Keep<-NA

for(x in 1:nrow(yesterday)){
	today$Current[today$Old==yesterday$Old[x]]<-yesterday$Current[x]
	today$Notes[today$Old==yesterday$Old[x]]<-yesterday$Notes[x]
	today$Keep[today$Old==yesterday$Old[x]]<-yesterday$Keep[x]
}

write.csv(today,file="Formicidae.combined.acc.nos.ABSENT.OLD.updated.csv")

stuff<-read.csv(file="~/Desktop/ants.ncbi.06july2016/Formicidae.combined.acc.nos.ABSENT.OLD.updated.csv",stringsAsFactors=FALSE,row.names=1,na.strings=c("","NA"))
stuff.keep<-stuff[stuff$Keep=="Yes",]
colnames(stuff.keep)<-c("Old","Taxon","Notes","Keep")
loci<-c("nuSSU","nuLSU","AbdA","COI","LR","Wg","EF1aF1","EF1aF2","ArgK","CAD","Top1","Ubx")
loci.tax<-c("Taxon","NCBI.ID",loci)
new.accs<-cbind(new[,loci.tax])
#pull taxa who's names are in antcat or their current names are missing from antcat from ncbi
things<-rbind(new.accs[new.accs$Taxon %in% old$Taxon,],new.accs[new.accs$Taxon %in% stuff.keep$Old,],new.accs[new.accs$Taxon %in% stuff.keep$Taxon,])

#then change names of those whose old names are found in current.missing.old (means they are missing from Antcat), and add their current name
for(p in 1:nrow(stuff.keep)){
	if(stuff.keep$Old[p]!=stuff.keep$Taxon[p]){
		things$NCBI.ID[things$Taxon==stuff.keep$Old[p]]<-NA
	}
	things$Taxon[things$Taxon==stuff.keep$Old[p]]<-stuff.keep$Taxon[p]
}
things$Taxon[duplicated(things$Taxon)]
gav<-things[!duplicated(things),]
gav$Taxon[duplicated(gav$Taxon)]

#now combine them

for(j in 1:nrow(gav)){
	if(is.na(gav$NCBI.ID[j])){
		for(locus in 1:length(loci)){
			if(isTRUE(is.na(gav[!is.na(gav$NCBI.ID) & gav$Taxon==gav$Taxon[j] & !is.na(gav[j,loci[locus]]),loci[locus]]))){
				print(loci[locus])
				print(gav[j,loci[locus]])
				gav[!is.na(gav$NCBI.ID) & gav$Taxon==gav$Taxon[j],loci[locus]]<-gav[j,loci[locus]]
			}
		}
	}
}

modified<-gav[!(is.na(gav$NCBI.ID) & duplicated(gav$Taxon)),]

old<-old[!duplicated(old),]

for(j in 1:nrow(modified)){
	if(modified$Taxon[j] %in% old$Taxon){
		modified$Species[j]<-old$Species[old$Taxon %in% modified$Taxon[j]]
		modified$Genus[j]<-old$Genus[old$Taxon %in% modified$Taxon[j]]
		modified$Tribe[j]<-old$Tribe[old$Taxon %in% modified$Taxon[j]]
		modified$Subfamily[j]<-old$Subfamily[old$Taxon %in% modified$Taxon[j]]
		modified$Palearctic[j]<-old$Palearctic[old$Taxon %in% modified$Taxon[j]]
		modified$Oceania[j]<-old$Oceania[old$Taxon %in% modified$Taxon[j]]
		modified$Neotropic[j]<-old$Neotropic[old$Taxon %in% modified$Taxon[j]]
		modified$Malagasy[j]<-old$Malagasy[old$Taxon %in% modified$Taxon[j]]
		modified$Australasia[j]<-old$Australasia[old$Taxon %in% modified$Taxon[j]]
		modified$Afrotropic[j]<-old$Afrotropic[old$Taxon %in% modified$Taxon[j]]
		modified$Indomalaya[j]<-old$Indomalaya[old$Taxon %in% modified$Taxon[j]]
	}
	if(!modified$Taxon[j] %in% old$Taxon){
		modified$Species[j]<-strsplit(modified$Taxon[j],"_")[[1]][2]
		modified$Genus[j]<-strsplit(modified$Taxon[j],"_")[[1]][1]
		modified$Tribe[j]<-old$Tribe[old$Genus %in% modified$Genus[j]][1]
		modified$Subfamily[j]<-old$Subfamily[old$Genus %in% modified$Genus[j]][1]
		modified$Palearctic[j]<-NA	
		modified$Oceania[j]<-NA
		modified$Neotropic[j]<-NA
		modified$Malagasy[j]<-NA
		modified$Australasia[j]<-NA
		modified$Afrotropic[j]<-NA
		modified$Indomalaya[j]<-NA
	}	
}

#list now updated with new names and geography
write.csv(modified,"ants.ready.to.roll.06july2016.csv",row.names=FALSE)

loci<-c("nuSSU","nuLSU","AbdA","COI","LR","Wg","EF1aF1","EF1aF2","ArgK","CAD","Top1","Ubx")
path<-"SPECIFYPATH"
filename<-"ants.ready.to.roll.06july2016.csv"
seq.refs<-read.csv(file=paste(path,filename,sep=""),stringsAsFactors=FALSE)
dir.create(paste(path,"/Formicidae",sep=""))

for(locus in loci){
	accs<-seq.refs[,locus][!is.na(seq.refs[,locus])]
	write(accs,file=paste(path,"Formicidae/",locus,".acc.nos.txt",sep=""))
}
