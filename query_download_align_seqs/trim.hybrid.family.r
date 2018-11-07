#This is called from within 6_blasting_ants_tax_hybrid_family2.py and uses R to parse the blast output

loci<-c("nuSSU", "nuLSU", "AbdA", "COI", "LR", "Wg", "EF1aF1", "EF1aF2", "ArgK", "CAD", "Top1", "Ubx")
clade<-"Formicidae"
folderpath<-"/FOLDERPATH/"
accfile<-"/FOLDERPATH/combined.updated.pruned.csv"
out.accfile<-"/FOLDERPATH/combined.updated.updated.pruned.csv"


orig.accfile<-read.csv(file=accfile,header=TRUE)
orig.accfile$Taxon<-gsub(" ","_",orig.accfile$Taxon)
write.csv(orig.accfile,file=out.accfile,row.names=FALSE)

for (locus in loci){
	bl<-read.table(file=as.character(paste(folderpath,clade,"/",locus,"/",locus,"_blast_for_edit_clean.txt",sep="")),sep="\t",header=FALSE,fill=TRUE,na.strings = "", quote="\"")
	colnames(bl)<-c("Query","LocalHit")
	bl$LocusHit<-bl$LocalHit
	#print(head(bl))
	for(l in loci){
		bl$LocusHit<-gsub(paste(l,".*",sep=""),l,bl$LocusHit)		
	}
	bl$LocalHit<-gsub(".*Formicidae.*","Formicidae",bl$LocalHit)
	bl$LocalHit<-gsub(".*Apidae.*","Outgroup",bl$LocalHit)
	bl$LocalHit<-gsub(".*Bethylidae.*","Outgroup",bl$LocalHit)
	bl$LocalHit<-gsub(".*Bradynobaenidae.*","Outgroup",bl$LocalHit)
	bl$LocalHit<-gsub(".*Crabronidae.*","Outgroup",bl$LocalHit)
	bl$LocalHit<-gsub(".*Mutillidae.*","Outgroup",bl$LocalHit)
	bl$LocalHit<-gsub(".*Pompilidae.*","Outgroup",bl$LocalHit)
	bl$LocalHit<-gsub(".*Sapygidae.*","Outgroup",bl$LocalHit)
	bl$LocalHit<-gsub(".*Scoliidae.*","Outgroup",bl$LocalHit)
	bl$LocalHit<-gsub(".*Sphecidae.*","Outgroup",bl$LocalHit)
	bl$LocalHit<-gsub(".*Tiphiidae.*","Outgroup",bl$LocalHit)
	bl$LocalHit<-gsub(".*Vespidae.*","Outgroup",bl$LocalHit)
	bl$LocalHit<-gsub(".*OUTGROUP.*","Outgroup",bl$LocalHit)
	bl$LocalHit<-gsub(".*Fungi.*","Fungi",bl$LocalHit)	#head(bl)
	write.csv(bl,file=as.character(paste(folderpath,clade,"/",locus,"/",locus,"_blast_for_edit_clean.csv",sep="")),row.names=FALSE)
	master<-read.csv(file=out.accfile,header=TRUE)
	master$Taxon<-as.character(master$Taxon)
	master$Taxon<-gsub(" ","_",master$Taxon)
	bl$Family<-NA
	for(p in 1:nrow(bl)){
		if(any(master$Taxon == bl$Query[p])){
			bl$Family[p]<-as.character(master$Family[master$Taxon == bl$Query[p]])
		}
	}
	pres.orig<-as.vector(master$Taxon[!is.na(master[locus])])
	cat(paste("\t...", length(pres.orig), "sequences in original", locus, "alignment...\n",sep=" "))
	good.new<-as.vector(bl$Query[bl$LocusHit==locus & bl$LocalHit==bl$Family])
	good.new<-good.new[!is.na(good.new)]
	cat(paste("\t\t...", length(good.new), "of", length(pres.orig), locus, "sequences passed BLAST...\n",sep=" "))
	in.old.not.new<-as.data.frame(setdiff(pres.orig,good.new))
	colnames(in.old.not.new)<-c("Taxon")
	if(nrow(in.old.not.new)>0){
		in.old.not.new$Accession<-NA
		in.old.not.new$LocalHit<-NA
		in.old.not.new$LocusHit<-NA
		for(j in 1:nrow(in.old.not.new)){
			if(in.old.not.new$Taxon[j] %in% bl$Query){
				in.old.not.new$LocalHit[j]<-bl$LocalHit[bl$Query %in% in.old.not.new$Taxon[j]]
				in.old.not.new$LocusHit[j]<-bl$LocusHit[bl$Query %in% in.old.not.new$Taxon[j]]
			}
			in.old.not.new$Accession[j]<-orig.accfile[,locus][orig.accfile$Taxon %in% in.old.not.new$Taxon[j]]
		}
	}
	cat(paste("\t\t...", nrow(in.old.not.new), "of", length(pres.orig), locus, "sequences FAILED BLAST...\n",sep=" "))
	master[locus][master$Taxon %in% in.old.not.new$Taxon,]<-NA
	master<-master[!rowSums(is.na(master[loci]))==length(loci),]
	write.csv(master,file=out.accfile,row.names=FALSE)
	write(good.new,file=as.character(paste(folderpath,clade,"/",locus,"/",locus,".taxa.blastok.txt",sep="")))
	#write(in.old.not.new,file=as.character(paste(folderpath,clade,"/",locus,"/",locus,".taxa.blastfailed.txt",sep="")))	
	write.csv(in.old.not.new,file=as.character(paste(folderpath,clade,"/",locus,"/",locus,".taxa.blastfailed.csv",sep="")),row.names=FALSE)
}
