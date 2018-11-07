#R script that takes tribe-level alignments and peforms profile alignments, ultimately resulting in a single alignment. Gblocks then used to trim alignment.

profile.align.tree.list.maker<-function(to.make,levels.for.tree,clade.level,locus,input.accfile){
	require(plyr)
	require(ape)
	input.accfile<-read.csv(input.accfile,stringsAsFactors=FALSE)
	mylist<-subset(input.accfile,select=c(levels.for.tree,locus))
	newlist<-mylist[!is.na(mylist[,locus]),]
	newlist<-as.data.frame(unclass(newlist))
	seqs.per.grp<-as.data.frame(table(newlist[clade.level]))
	colnames(seqs.per.grp)<-c(as.character(clade.level),"Abund")
	newlist<-join(newlist,seqs.per.grp,by=as.character(clade.level))
	newlist.b<-newlist
	newlist.b[,locus]<-NULL
	tree.dat<-unique(newlist.b)
	tax.hier<-paste("~",gsub(" ","/",paste(levels.for.tree,collapse=" ")),sep="")
	tr<-as.phylo(as.formula(tax.hier),data=tree.dat)
	tr.ord<-reorder(tr,"postorder")
	org<-as.data.frame(tr.ord$edge)
	colnames(org)<-c("From","To")
	nn<-unique(org$From)
	nn.fr<-plyr::count(org$From)
	colnames(nn.fr)<-c("Node.No","Desc")
	nn.fr<-nn.fr[match(nn.fr$Node.No,nn),]
	org$To.Name<-tr.ord$tip.label[tr.ord$edge[,2]]
	org$To.Name[is.na(org$To.Name)]<-org$To[is.na(org$To.Name)]
	for(tax in 1:nrow(org)){
		org$Abund[tax]<-newlist$Abund[newlist[clade.level]==org$To.Name[tax]][1]
	}
	if(to.make=="tree"){
		return(tr.ord)
	}	
	if(to.make=="large.summary"){
		return(org)
	}	
	if(to.make=="par.desc"){
		return(nn.fr)
	}
}


profile.aligner.merge<-function(clade,clade.level,locus,organization,node.nos.fr,folderpath,path.to.mafft,path.to.macse,path.to.gblocks,dat.type,threads.mafft,amb.rem,translation.cat,tax.stand=NULL){
	require(ape)
	require(seqinr)
	for(z in 1:nrow(node.nos.fr)){
		nds<-organization$To.Name[organization$From==node.nos.fr$Node.No[z]]
		#cat.files<-as.vector(paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",nds,"_",locus,"_aligned.fasta",sep=""))
		cat.files<-as.character(paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",nds,"_",locus,"_aligned.fasta",sep="", collapse=" "))
		cat.outfile<-as.character(paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",node.nos.fr$Node.No[z],"_",locus,"_cat_aligned.fasta",sep=""))
		cat.commands<-as.character(paste("cat",paste(cat.files,sep=" "),">",cat.outfile,sep=" "))
		system(cat.commands,intern=TRUE)
		ruby.table.maker.file<-as.character(paste(folderpath,"makemergetable.rb",sep=""))
		sub.msa.table<-as.character(paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",node.nos.fr$Node.No[z],"_",locus,"_submsa_table",sep=""))
		ruby.commands<-as.character(paste("ruby",ruby.table.maker.file, cat.files, ">",sub.msa.table,sep=" "))
		system(ruby.commands,intern=TRUE)
		output.file<-as.character(paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",node.nos.fr$Node.No[z],"_",locus,"_aligned.fasta",sep=""))	
		#mafft.commands<-paste(path.to.mafft,"--localpair --maxiterate 1000 --preservecase --merge", sub.msa.table, cat.outfile, ">",output.file, sep=" ")
		#make something here to get number of seqs going into alignment
		#print(cat.files)
		#this is problematic if both files are joined together like this "file1 file2" instead of "file1" "file2"
		#ntax<-sum(sapply(cat.files,function(zaz) length(read.fasta(file=zaz,seqtype=dat.type[ind.loc],forceDNAtolower=FALSE))))
		#so instead count number in cat.outfile
		if(dat.type!="AA"){
			ntax<-length(read.fasta(file=cat.outfile,seqtype="DNA",forceDNAtolower=FALSE))	
		}
		if(dat.type=="AA"){
			ntax<-length(read.fasta(file=cat.outfile,seqtype="AA",forceDNAtolower=FALSE))				
		}
		#then make small decision tree for mafft commands - if less than 200 seqs, use G-INS-i (and --leavegappyregion), but if over, use FFT-NS-i (w/o --leavegappyregion)
		#--leavegappyregion only works well in small alignments (Katoh & Standley 2016), so remove from large aligns
		#does not seem wise to set unalignlevel (VSM) higher than 0.8 (Katoh & Standley 2016).
		if(ntax<200){
			#mafft.commands<-paste(path.to.mafft,"--globalpair --maxiterate 1000 --thread", threads, "--preservecase --merge", sub.msa.table, cat.outfile, ">",output.file, sep=" ")
			#did G-INS-i instead
			#mafft.commands<-paste(path.to.mafft,"--globalpair --maxiterate 1000 --thread", threads, "--preservecase --unalignlevel 0.8 --leavegappyregion --merge", sub.msa.table, cat.outfile, ">", output.file, sep=" ")	
			#did E-INS-i instead
			mafft.commands<-paste(path.to.mafft,"--genafpair --maxiterate 1000 --thread", threads.mafft, "--preservecase --merge", sub.msa.table, cat.outfile, ">", output.file, sep=" ")	
		}
		if(ntax>=200){
			#did FFT-NS-i...cannot add unalignlevel as this is only available when using G-INS-i
			mafft.commands<-paste(path.to.mafft,"--retree 2 --maxiterate 1000 --thread", threads.mafft, "--preservecase --merge", sub.msa.table, cat.outfile, ">",output.file, sep=" ")	
		}
		system(mafft.commands,intern=FALSE)
		Sys.sleep(1)
		#Refine alignment w MACSE
		if(dat.type=="codon"){
			output.file<-as.character(paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",node.nos.fr$Node.No[z],"_",locus,"_aligned.fasta",sep=""))	
			output.file.unrefined<-as.character(paste(strsplit(output.file,split=".fasta")[[1]],"_unrefined.fasta",sep=""))
			copy.commands.first<-paste("cp", output.file, output.file.unrefined, sep=" ")
			system(copy.commands.first)
			commands.to.macse<-paste("java -jar -Xmx1200m", path.to.macse,"-prog refineAlignment -align", output.file.unrefined, "-gc_def", translation.cat, "-out_AA", paste(output.file,"_AA.fasta",sep=""), "-out_NT", paste(output.file,"_NT.fasta",sep=""), "-optim 1 -ext_gap_ratio 0.0001 -gap_op 1", sep=" ")
			system(commands.to.macse)
			copy.commands<-paste("cp", as.character(paste(output.file,"_NT.fasta",sep="")), output.file, sep=" ")
			system(copy.commands)
			commands.to.perl<-paste("perl -i -pe 's/[!]/-/g'", output.file, sep=" ")
			system(commands.to.perl)
			Sys.sleep(1)
		}	
	}
	rename.commands<-paste("cp", paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",node.nos.fr$Node.No[nrow(node.nos.fr)],"_",locus,"_aligned.fasta",sep=""), paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""), sep=" ")
	system(rename.commands,intern=TRUE)
	Sys.sleep(1)
	#Remove taxon standards if included previously
	if(!is.null(tax.stand)){	
		if(dat.type!="AA"){
			mat<-read.fasta(file=paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""),seqtype="DNA",forceDNAtolower=FALSE)
			mat<-mat[-grep("Standard_.*",names(mat))]
			write.fasta(sequences=mat,names=names(mat),file=paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""),nbchar=1000000)
		}
		if(dat.type=="AA"){
			mat<-read.fasta(file=paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""),seqtype="AA",forceDNAtolower=FALSE)		
			mat<-mat[-grep("Standard_.*",names(mat))]
			write.fasta(sequences=mat,names=names(mat),file=paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""),nbchar=1000000)
		}
	}
	#added condition that gblocks only used if alignment specified
	if(amb.rem=="y"){
		if(dat.type!="AA"){
			full.align<-read.fasta(file=paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""),seqtype="DNA",forceDNAtolower=FALSE)
		}
		if(dat.type=="AA"){
			full.align<-read.fasta(file=paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""),seqtype="AA",forceDNAtolower=FALSE)
		}
		ntax<-length(full.align)
		if(dat.type=="DNA"){
			t.var<-as.character(paste("-t=","d",sep=""))
		}
		if(dat.type=="AA"){
			t.var<-as.character(paste("-t=","p",sep=""))
		}
		if(dat.type=="codon"){
			t.var<-as.character(paste("-t=","c",sep=""))
		}
		subdat<-ceiling((ntax/2)+0.5)
		b2.var<-as.character(paste("-b2=",subdat,sep=""))
		b4.var<-as.character(paste("-b4=","5",sep=""))	
		b5.var<-as.character(paste("-b5=","h",sep=""))
		file.extension<-as.character(paste("-e=","-gb5",sep=""))	
		gblocks.commands<-paste(path.to.gblocks, paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""), t.var, b2.var, b4.var, b5.var, file.extension, sep=" ")
		print(gblocks.commands)
		system(gblocks.commands,intern=TRUE)
		Sys.sleep(1)
		b4.var<-as.character(paste("-b4=","2",sep=""))
		file.extension<-as.character(paste("-e=","-gb2",sep=""))	
		gblocks.commands<-paste(path.to.gblocks, paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta",sep=""), t.var, b2.var, b4.var, b5.var, file.extension, sep=" ")
		cat(gblocks.commands)
		system(gblocks.commands,intern=TRUE)
		Sys.sleep(1)
		#to.rename<-read.fasta(file=as.character(paste(output.file,"-gb",sep="")),seqtype=dat.type[ind.loc],forceDNAtolower=FALSE)
		#write.fasta(sequences=to.rename,names=names(to.rename),file=output.file,nbchar=1000000)
		#rename.commands<-paste("cp", paste(output.file,"-gb",sep=""), output.file, sep=" ")
		#system(rename.commands,intern=TRUE)
		#Sys.sleep(1)
	}
	if(dat.type!="AA"){
			full.align<-read.fasta(file=paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta-gb5",sep=""),seqtype="DNA",forceDNAtolower=FALSE)
	}
	if(dat.type=="AA"){
			full.align<-read.fasta(file=paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta-gb5",sep=""),seqtype="AA",forceDNAtolower=FALSE)
	}
	outs<-potential.outgroups[potential.outgroups %in% names(full.align)]
	outs<-paste(outs,collapse=",")
	print(outs)
	#raxml.commands<-paste(path.to.raxml, "-T", threads, "-f a -x 12345 -p 12345 -o", outs, "-m GTRCAT -# 1000 -s", paste(folderpath,clade,"/",locus,"/",clade.level,"_Matrices","/",locus,"_aligned.fasta-gb",sep=""),"-n", paste(locus,"_RAxML",sep=""))
	#print(raxml.commands)
	#system(raxml.commands,intern=TRUE)
	#Sys.sleep(1)
}

profile.to.end.single<-function(levels.for.tree=tree.levs,locus=locus,input.accfile=input.accfile,tax.stands=NULL,clade=clade,clade.level=clade.level,organization=ind.org,node.nos.fr=ind.node.no.fr,folderpath=folderpath,path.to.mafft=path.to.mafft, path.to.macse=path.to.macse,threads.mafft=threads.mafft,amb.rem="y",dat.type=dat.type,translation.cat=translation.cat,path.to.gblocks=path.to.gblocks,seq.refs=seq.refs,potential.outgroups=potential.outgroups){
	ind.org<-profile.align.tree.list.maker(to.make="large.summary",levels.for.tree=tree.levs,clade.level=clade.level,locus=locus,input.accfile=input.accfile)
	ind.node.no.fr<-profile.align.tree.list.maker(to.make="par.desc",levels.for.tree=tree.levs,clade.level=clade.level,locus=locus,input.accfile=input.accfile)
	profile.aligner.merge(clade=clade,clade.level=clade.level,locus=locus,dat.type=dat.type,translation.cat=translation.cat,organization=ind.org,node.nos.fr=ind.node.no.fr,folderpath=folderpath,path.to.mafft=path.to.mafft,path.to.macse=path.to.macse,path.to.gblocks=path.to.gblocks,threads.mafft=threads.mafft,amb.rem="y")
}

loci<-c("nuSSU","nuLSU","AbdA","COI","LR","Wg","EF1aF1","EF1aF2","ArgK","CAD","Top1","Ubx")
dat.type<-c("DNA","DNA","codon","codon","codon","codon","codon","codon","codon","codon","codon","codon")
translation.cat<-c(0,0,1,5,1,1,1,1,1,1,1,1)
folderpath<-"/PATHTOFILES/"
clade="Hymenoptera"
clade.level<-"Tribe"
path.to.mafft<-"/PATHTOMAFFT/mafft"
path.to.macse<-"/PATHTOMACSE/macse_v1.2.jar"
tree.levs<-c("Family","Subfamily","Tribe")
path.to.gblocks<-"/PATHTOGBLOCKS/Gblocks"
threads<-20
path.to.raxml<-"/PATHTORAXML/raxmlHPC-PTHREADS-SSE3"
input.accfile<-as.character(paste(folderpath,clade,"/","combined.updated.updated.pruned.wouts.csv",sep=""))
seq.refs<-read.csv(file=as.character(input.accfile),stringsAsFactors=FALSE)
potential.outgroups<-c(seq.refs$Taxon[seq.refs$Family=="Bethylidae"])
taxon.standards<-c("Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera","Apis_mellifera")

if(any(dat.type!="codon")){
	dat.type.non.codon<-dat.type[dat.type!="codon"]
	loci.non.codon<-loci[dat.type!="codon"]
	translation.cat.non.codon<-translation.cat[dat.type!="codon"]
	if(!missing(taxon.standards)){
		taxon.standards.codon<-taxon.standards[dat.type!="codon"]
	}
	if(missing(taxon.standards)){
		taxon.standards.codon<-NULL		
	}
	for(ind.loc in 1:length(loci.non.codon)){
		profile.to.end.single(levels.for.tree=tree.levs,locus=loci.non.codon[ind.loc],taxon.standards.codon[ind.loc],input.accfile=input.accfile,clade=clade,clade.level=clade.level,organization=ind.org,node.nos.fr=ind.node.no.fr,folderpath=folderpath,path.to.mafft=path.to.mafft, threads.mafft=threads,amb.rem="y",dat.type=dat.type.non.codon[ind.loc],translation.cat=translation.cat.non.codon[ind.loc],path.to.gblocks=path.to.gblocks,seq.refs=seq.refs,potential.outgroups=potential.outgroups)
	}
}

if(any(dat.type=="codon")){
	require(doParallel)
	require(seqinr)
	dat.type.codon<-dat.type[dat.type=="codon"]
	loci.codon<-loci[dat.type=="codon"]
	translation.cat.codon<-translation.cat[dat.type=="codon"]
	if(!missing(taxon.standards)){
		taxon.standards.codon<-taxon.standards[dat.type=="codon"]
	}
	if(missing(taxon.standards)){
		taxon.standards.codon<-NULL		
	}
	cl<-makeCluster(20)
	length(cl)
	registerDoParallel(cl)
	foreach(qv=1:length(loci.codon)) %dopar% profile.to.end.single(levels.for.tree=tree.levs,locus=loci.codon[qv],tax.stands=taxon.standards.codon[qv],input.accfile=input.accfile,clade=clade,clade.level=clade.level,organization=ind.org,node.nos.fr=ind.node.no.fr,folderpath=folderpath,path.to.mafft=path.to.mafft, threads.mafft=1,amb.rem="y",dat.type=dat.type.codon[qv],translation.cat=translation.cat.codon[qv],path.to.gblocks=path.to.gblocks,path.to.macse=path.to.macse,seq.refs=seq.refs,potential.outgroups=potential.outgroups)
}
stopCluster(cl)
