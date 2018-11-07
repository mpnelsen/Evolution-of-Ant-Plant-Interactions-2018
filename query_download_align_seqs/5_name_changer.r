#This R file reads in alignments and updates names of taxa to be consistent with those in spreadsheet

#be sure to first copy and rename 5.8s output file ITS_Parsed.fasta

require(seqinr)
loci<-c("nuSSU","nuLSU","AbdA","COI","LR","Wg","EF1aF1","EF1aF2","ArgK","CAD","Top1","Ubx")
dat.type<-c("DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA","DNA")
seq.refs<-read.csv(file="/PATHTO/ants.ready.to.roll.06july2016.csv",stringsAsFactors=FALSE)
path<-"/PATH"

for(locus in loci[3:12]){
	align<-read.fasta(file=paste(path,"/",locus,"/",locus,"_CDS_Parsed.fasta",sep=""),seqtype="DNA",forceDNAtolower=FALSE)
	if(locus=="COI"){
		g.align<-read.fasta(file=paste(path,"/",locus,"/",locus,"_GENE_Parsed.fasta",sep=""),seqtype="DNA",forceDNAtolower=FALSE)
		combo.align<-c(g.align,align)
		align<-combo.align
	}
	write.fasta(sequences=align,names=names(align),file.out=paste(path,"/",locus,"/",locus,"_Parsed.fasta",sep=""),nbchar=1000000)
}


for(locus in loci){
	require(ape)
	new.align<-read.fasta(file=paste(path,"/",locus,"/",locus,"_Parsed.fasta",sep=""),seqtype=dat.type[loci==locus],forceDNAtolower=FALSE)
	names(new.align)<-gsub("[.].*","",names(new.align))
	for(p in 1:length(names(new.align))){
		#print(seq.refs$Taxon[seq.refs[,locus] %in% names(new.align)[p]])
		#print(length(new.align))
		#print(names(new.align)[duplicated(names(new.align))])
		#print(names(new.align)[p])
		names(new.align)[p]<-seq.refs$Taxon[seq.refs[,locus] %in% names(new.align)[p]]
	}
	write.fasta(sequences=new.align,names=gsub(" ","_",names(new.align)),file=paste(path,"/",locus,"/",locus,"_Renamed_Parsed.fasta",sep=""),nbchar=1000000)
	names.w.accs<-seq.refs$Taxon[!is.na(seq.refs[,locus])]
	names.missing<-names.w.accs[!names.w.accs %in% names(new.align)]
	print(locus)
	print(names.missing)
	#print(seq.refs[,locus][names.missing %in% seq.refs$Taxon])
	seq.refs[seq.refs$Taxon %in% names.missing,locus]<-NA
	#line below will remove taxa lacking any ribosomal sequences
	#seq.refs<-seq.refs[rowSums(is.na(seq.refs[,loci]))!=length(loci),]
}

seq.refs$Taxon[is.na(seq.refs$nuSSU) & is.na(seq.refs$nuLSU) & is.na(seq.refs$AbdA) & is.na(seq.refs$Wg) & is.na(seq.refs$LR) & is.na(seq.refs$COI) & is.na(seq.refs$EF1aF1) & is.na(seq.refs$EF1aF2) & is.na(seq.refs$ArgK) & is.na(seq.refs$CAD) & is.na(seq.refs$Top1) & is.na(seq.refs$Ubx)]
write.csv(seq.refs,file=paste(path,"/","combined.updated.pruned.csv",sep=""),row.names=FALSE)


#####MANUALLY ADD IN FAMILY COLUMN TO COMBINED.UPDATED.PRUNED.CSV#####
