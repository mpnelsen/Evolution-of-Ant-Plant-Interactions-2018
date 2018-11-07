#R function/wrapper to align sequences

initial.aligner<-function(locus,seq.type,trans.cat,folderpath,threads,input.accfile,path.to.mafft,path.to.macse=NULL,tax.acc.list=NULL,clade,clade.level,tax.stand=NULL){
	require(seqinr)
	require(doParallel)
	#tax.acc.list[clade.level][is.na(tax.acc.list[clade.level])]<-"IncertaeSedis"
	dir.create(as.character(paste(folderpath,clade,"/",locus,sep="")))
	dir.create(as.character(paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices",sep="")))
	if(seq.type=="AA"){
		big.matrix<-read.fasta(file=paste(folderpath,clade,"/",locus,"_Renamed_wouts_Parsed.fasta",sep=""),seqtype="AA",forceDNAtolower=FALSE)
	}
	if(seq.type!="AA"){
		big.matrix<-read.fasta(file=paste(folderpath,clade,"/",locus,"_Renamed_wouts_Parsed.fasta",sep=""),seqtype="DNA",forceDNAtolower=FALSE)
	}
	if(any(seq.type==c("DNA","AA"))){	
		clades.unique<-unique(tax.acc.list[clade.level][!is.na(tax.acc.list[locus])])
		cl.no<-length(clades.unique)
		for(sub.clade in clades.unique){
			reduced.matrix<-big.matrix[names(big.matrix) %in% tax.acc.list$Taxon[!is.na(tax.acc.list[locus]) & tax.acc.list[clade.level]==sub.clade]]
			#if standards included, re-select w standard and rename it.
			if(!is.null(tax.stand)){
				print(tax.stand)
				reduced.matrix<-big.matrix[names(big.matrix) %in% tax.acc.list$Taxon[!is.na(tax.acc.list[locus]) & tax.acc.list[clade.level]==sub.clade] | names(big.matrix) %in% tax.stand]
				print(length(reduced.matrix))
				#change name of standard to be Standard_sub.clade_tax.stand
				names(reduced.matrix)[names(reduced.matrix) %in% tax.stand]<-paste("Standard",sub.clade,tax.stand,sep="_")
				print(names(reduced.matrix))
			}
			if(length(names(reduced.matrix))==1){
				write.fasta(sequences=reduced.matrix,names=names(reduced.matrix),file=paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",sub.clade,"_",locus,"_aligned.fasta",sep=""),nbchar=1000000)
			}
			if(length(names(reduced.matrix))>1){
				cat(paste("\t...Aligning", sub.clade, locus, "sequences...\n",sep=" "))
				write.fasta(sequences=reduced.matrix,names=names(reduced.matrix),file=paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",sub.clade,"_",locus,"_unaligned.fasta",sep=""),nbchar=1000000)
				input.file<-as.character(paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",sub.clade,"_",locus,"_unaligned.fasta",sep=""))
				output.file<-as.character(paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",sub.clade,"_",locus,"_aligned.fasta",sep=""))	
				if(length(names(reduced.matrix))<200){
					#do G-INS-i
					#mafft.commands<-paste(path.to.mafft,"--globalpair --maxiterate 1000 --thread", threads, "--preservecase --unalignlevel 0.8 --leavegappyregion", input.file, ">", output.file, sep=" ")	
					#do E-INS-i
					mafft.commands<-paste(path.to.mafft,"--genafpair --maxiterate 1000 --thread", threads, "--preservecase", input.file, ">", output.file, sep=" ")					
				}
				if(length(names(reduced.matrix))>=200){
					#did FFT-NS-i...cannot add unalignlevel as this is only available when using G-INS-i
					mafft.commands<-paste(path.to.mafft,"--retree 2 --maxiterate 1000 --thread", threads, "--preservecase", input.file, ">",output.file, sep=" ")	
				}
				print(mafft.commands)
				system(mafft.commands,intern=TRUE)
				Sys.sleep(0.25)
			}
		}
	}		
	if(seq.type=="codon"){
		clades.unique<-unique(tax.acc.list[clade.level][!is.na(tax.acc.list[locus])])
		cl.no<-length(clades.unique)
		cl<-makeCluster(threads)
		registerDoParallel(cl)
		macse.aligner<-function(bm,sc,loc,cld,cldlev,fp,tal,tc,ts=NULL){
			reduced.matrix<-bm[names(bm) %in% tal$Taxon[!is.na(tal[loc]) & tal[cldlev]==sc]]
			if(!is.null(ts)){
				print(ts)
				reduced.matrix<-bm[names(bm) %in% tal$Taxon[!is.na(tal[loc]) & tal[cldlev]==sc] | names(bm) %in% ts]
				print(length(reduced.matrix))
				#change name of standard to be Standard_sub.clade_tax.stand
				names(reduced.matrix)[names(reduced.matrix) %in% ts]<-paste("Standard",sc,ts,sep="_")
				print(names(reduced.matrix))
			}
			seqinr::write.fasta(sequences=reduced.matrix,names=names(reduced.matrix),file=paste(fp,cld,"/",loc,"/",cldlev,"_Matrices","/",sc,"_",loc,"_unaligned.fasta",sep=""),nbchar=1000000)
			cat(paste("\t...Aligning", sc, loc, "sequences...\n",sep=" "))
			input.file<-as.character(paste(fp,cld,"/",loc,"/",cldlev,"_Matrices","/",sc,"_",loc,"_unaligned.fasta",sep=""))
			output.file<-as.character(paste(fp,cld,"/",loc,"/",cldlev,"_Matrices","/",sc,"_",loc,"_aligned.fasta",sep=""))	
			commands.to.macse<-paste("java -jar -Xmx1200m", path.to.macse,"-prog alignSequences -ext_gap_ratio 0.0001 -gap_op 1 -seq", input.file, "-gc_def", tc, "-out_AA", paste(output.file,"_unrefined_AA.fasta",sep=""), "-out_NT", paste(output.file,"_unrefined_NT.fasta",sep=""), sep=" ")
			system(commands.to.macse)
			print(commands.to.macse)
			#Decided not to refine
			#Refine alignment w MACSE
			#cat(paste("\t\t...Refining", sc, loc, "alignment...\n",sep=" "))
			#commands.to.macse<-paste("java -jar -Xmx1200m", path.to.macse,"-prog refineAlignment -align", paste(output.file,"_unrefined_NT.fasta",sep=""), "-gc_def", tc, "-out_AA", paste(output.file,"_refined_AA.fasta",sep=""), "-out_NT", paste(output.file,"_refined_NT.fasta",sep=""), "-optim 1 -ext_gap_ratio 0.0001 -gap_op 1", sep=" ")
			#system(commands.to.macse)
			copy.commands<-paste("cp", as.character(paste(output.file,"_unrefined_NT.fasta",sep="")), output.file, sep=" ")
			system(copy.commands)
			commands.to.perl<-paste("perl -i -pe 's/[!]/-/g'", output.file, sep=" ")
			system(commands.to.perl)	
		}	
		foreach(x=clades.unique) %dopar% macse.aligner(bm=big.matrix,sc=x,loc=locus,cld=clade,cldlev=clade.level,fp=folderpath,tal=tax.acc.list,tc=trans.cat,ts=tax.stand)
		stopCluster(cl)
	}
}



clade="Hymenoptera"
clade.level<-"Tribe"
loci<-c("nuSSU","nuLSU","AbdA","COI","LR","Wg","EF1aF1","EF1aF2","ArgK","CAD","Top1","Ubx")
dat.type<-c("DNA","DNA","codon","codon","codon","codon","codon","codon","codon","codon","codon","codon")
translation.cat<-c(0,0,1,5,1,1,1,1,1,1,1,1)
folderpath<-"/PATHFORFILES/"
threads<-20
input.accfile<-as.character(paste(folderpath,clade,"/","combined.updated.updated.pruned.wouts.csv",sep=""))
path.to.mafft<-"/PATHTOMAFFT/mafft"
path.to.macse<-"/PATHTOMACSE/macse_v1.2.jar"
tax.acc.list<-read.csv(file=input.accfile,stringsAsFactors=FALSE)
taxon.standards<-c("Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera")

for(g in 1:length(loci)){
	initial.aligner(locus=loci[g],seq.type=dat.type[g],trans.cat=translation.cat[g],tax.stand=taxon.standards[g],folderpath=folderpath,threads=threads,input.accfile=input.accfile,path.to.mafft=path.to.mafft,path.to.macse=path.to.macse,tax.acc.list=tax.acc.list,clade=clade,clade.level=clade.level)
}

