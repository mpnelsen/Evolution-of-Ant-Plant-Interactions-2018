library(traits)
library(taxize)
library(rentrez)
library(plyr)

grab.acc.nos.from.clade<-function(clade.to.search,locus,search.terms,excludes=NULL){
	require(traits)
	require(taxize)
	require(rentrez)
	require(plyr)
	#require(httr)
	#httr::set_config(timeout(seconds = 30))
	#httr::timeout(seconds = 30)
	#found nice trick to retry if fail here:
	#http://stackoverflow.com/questions/20770497/how-to-retry-a-statement-on-error
	search.res<-NULL
	attempt<-1
	while(is.null(search.res) && attempt<=5){
		if(attempt>1){
			Sys.sleep(runif(1,min=1,max=10))
		}
		cat(paste("\t...Searching for", clade.to.search, locus, "sequences...Attempt", attempt, "\n",sep=" "))
		try(search.res<-ncbi_searcher(taxa=clade.to.search,seqrange="350:5000",limit=40000,verbose=FALSE,fuzzy=TRUE,entrez_query=paste(search.terms,excludes,sep=" ")))
		attempt<-attempt+1
	}
	#search.res<-ncbi_searcher(taxa=clade.to.search,seqrange="350:5000",limit=50000,verbose=FALSE,fuzzy=TRUE,entrez_query=search.terms)
	#print(search.res[,"taxon"])
	cat(paste("\t\t...Found", nrow(search.res), locus, "sequences...\n", sep=" "))
	#print(search.res[,c("taxon","acc_no")])
	if(nrow(search.res)>=50000){
		cat(paste("\t\t\t...", locus, " SEARCH LIMIT REACHED (50,000 accessions) - SEARCH HALTED = try searching a smaller clade...\n", sep=""))
		search.res.clean<-data.frame(NA,NA)
		colnames(search.res.clean)<-c("NCBI.ID",locus)
		return(search.res.clean)			
	}
	if(nrow(search.res)==0){
		cat(paste("\t\t...No",clade.to.search,locus,"accessions found...\n"))		
		search.res.clean<-data.frame(NA,NA)
		colnames(search.res.clean)<-c("NCBI.ID",locus)
		return(search.res.clean)			
	}
	#function has trouble parsing acc no properly if NG_ seqs are retrieved -> strips the NG_ and left w numbers - see NG_013161...but GI is fine
	#which will only list numbers
	#these will grab weird things w genbak2uid
	#so subset and make sure things like this are removed
	else{
		#This removes if there is no letters...so retain to keep reference seqs
		search.res<-search.res[grepl("[A-Z]",search.res$acc_no),]
		cat(paste("\t\t\t...Filtering", clade.to.search, "environmental and unidentified...\n",sep=" "))
		search.res<-search.res[!grepl("[0-9]|sp[.]|aff[.]|Cf[.]|cf[.]| cf |var[.]|subsp[.]|f[.]|auct[.]|UNVERIFIED",search.res$taxon),]
		if(nrow(search.res)==0){
			cat(paste("\t\t\t...No",clade.to.search,locus,"accessions found after filtering...\n"))		
			search.res.clean<-data.frame(NA,NA)
			colnames(search.res.clean)<-c("NCBI.ID",locus)
			return(search.res.clean)			
		}	
		#print(search.res[,c("taxon","acc_no")])
		dups.first<-unique(search.res$taxon[duplicated(search.res$taxon)])
		cat(paste("\t\t\t...Filtering", clade.to.search, "duplicate names...\n",sep=" "))
		if(length(dups.first)>0){
			dup.vals.first<-unique(search.res$taxon[duplicated(search.res$taxon)])
			#print(dup.vals.first)
			for(q in 1:length(dup.vals.first)){
				#add line here to preferentially retain reference seqs
				if(any(search.res$taxon==dup.vals.first[q] & grepl("^[A-Z][A-Z]_",search.res$acc_no))){
					search.res[search.res$taxon==dup.vals.first[q] & !grepl("^[A-Z][A-Z]_",search.res$acc_no),]<-NA
					search.res<-search.res[!is.na(search.res$taxon),]
					#print(search.res)
				}		
			}
		}
		#print(search.res[,c("taxon","acc_no")])
		cat(paste("\t\t...Retained", nrow(search.res), "taxa with", locus, "sequences after removing duplicate names...\n", sep=" "))
		dups.second<-unique(search.res$taxon[duplicated(search.res$taxon)])
		cat(paste("\t\t...Checking", clade.to.search, "for refseqs...\n",sep=" "))
		if(length(dups.second)>0){
			dup.vals.second<-unique(search.res$taxon[duplicated(search.res$taxon)])
			for(q in 1:length(dup.vals.second)){
				max.length.for.dup.second<-search.res$acc_no[search.res$taxon==dup.vals.second[q] & search.res$length==max(search.res$length[search.res$taxon==dup.vals.second[q]])][1]
				search.res[search.res$taxon==dup.vals.second[q] & search.res$acc_no!=max.length.for.dup.second,]<-NA
				search.res<-search.res[!is.na(search.res$taxon),]			
			}
		}
		#print(search.res[,c("taxon","acc_no")])
		search.res$NCBI.ID<-NA
		cat(paste("\t\t...Fetching taxid numbers for", clade.to.search, locus, "accessions...\n",sep=" "))
		#need to break this up to avoid timeouts - break into groups of 500
		#no.hits<-nrow(search.res)
		if(nrow(search.res)<=500){
			attempt<-1
			while(is.na(search.res$NCBI.ID[1]) && attempt<=5){
				if(attempt>1){
					Sys.sleep(runif(1,min=1,max=10))
				}
				try(search.res$NCBI.ID<-genbank2uid(search.res$acc_no[1:nrow(search.res)]))
				attempt<-attempt+1
			}		
		} else {
			for(p in 1:ceiling(nrow(search.res)/500)){
				if(p<ceiling(nrow(search.res)/500)){
					attempt<-1
					while(is.na(search.res$NCBI.ID[((500*p)-499)]) && attempt<=5){
						if(attempt>1){
							Sys.sleep(runif(1,min=1,max=10))
						}	
						try(search.res$NCBI.ID[((500*p)-499):(500*p)]<-genbank2uid(search.res$acc_no[((500*p)-499):(500*p)]))
						attempt<-attempt+1
					}	
				} else {
					attempt<-1
					while(is.na(search.res$NCBI.ID[((500*p)-499)]) && attempt<=5){
						if(attempt>1){
							Sys.sleep(runif(1,min=1,max=10))
						}
						try(search.res$NCBI.ID[((500*p)-499):nrow(search.res)]<-genbank2uid(search.res$acc_no[((500*p)-499):nrow(search.res)]))
						attempt<-attempt+1
					}
				}
			}
		}
		dups.third<-unique(search.res$NCBI.ID[duplicated(search.res$NCBI.ID)])
		#print(dups.third)
		cat(paste("\t\t...Checking", clade.to.search, "taxids for refseqs again...\n",sep=" "))
		if(length(dups.third)>0){
			dup.vals.third<-unique(search.res$NCBI.ID[duplicated(search.res$NCBI.ID)])
			#print(dup.vals.third)
			#print(head(search.res))
			for(q in 1:length(dup.vals.third)){
				if(any(search.res$NCBI.ID==dup.vals.third[q] & grepl("^[A-Z][A-Z]_",search.res$acc_no))){
					search.res[search.res$NCBI.ID==dup.vals.third[q] & !grepl("^[A-Z][A-Z]_",search.res$acc_no),]<-NA
					search.res<-search.res[!is.na(search.res$NCBI.ID),]			
				}
			}
		}
		#print(search.res[,c("taxon","acc_no")])
		dups.fourth<-unique(search.res$NCBI.ID[duplicated(search.res$NCBI.ID)])
		cat(paste("\t\t\t...Filtering", clade.to.search, "duplicate taxids...\n",sep=" "))
		if(length(dups.fourth)>0){
			dup.vals.fourth<-unique(search.res$NCBI.ID[duplicated(search.res$NCBI.ID)])
			for(q in 1:length(dup.vals.fourth)){
				max.length.for.dup.fourth<-search.res$acc_no[search.res$NCBI.ID==dup.vals.fourth[q] & search.res$length==max(search.res$length[search.res$NCBI.ID==dup.vals.fourth[q]])][1]
				search.res[search.res$NCBI.ID==dup.vals.fourth[q] & search.res$acc_no!=max.length.for.dup.fourth,]<-as.character(NA)
				search.res<-search.res[!is.na(search.res$NCBI.ID),]			
			}
		}
	if(nrow(search.res)==0){
		cat(paste("\t\t...No",clade.to.search,locus,"accessions found after filtering steps...\n"))		
		search.res.clean<-data.frame(NA,NA)
		colnames(search.res.clean)<-c("NCBI.ID",locus)
		return(search.res.clean)			
	}
	search.res.clean<-subset(search.res, select=c(NCBI.ID,acc_no))
	colnames(search.res.clean)<-c("NCBI.ID",locus)
	cat(paste("\t\t...Retained", nrow(search.res.clean), "taxa with", locus, "sequences after removing duplicate taxids...\n", sep=" "))
	return(search.res.clean)	
	}
}

get.accs.for.genes.in.clade<-function(clade.to.search,locus.list,search.terms){
	require(plyr)
	cat(paste("...Searching NCBI for", clade.to.search, "...\n",sep=" "))
	for(ll in 1:length(locus.list)){
		res.search<-grab.acc.nos.from.clade(clade.to.search=clade.to.search,locus=locus.list[ll],search.terms=search.terms[ll],excludes=excludes)
		if(ll==1){
			all.search<-res.search
		}
		if(ll>1){
			all.search<-join(all.search,res.search,by="NCBI.ID",type="full")
		}
	}
	all.search$Taxon<-NA
	#this will then get rid of any that totally lack an NCBI.ID, which might happen if your search for a gene returned nothing at all
	all.search<-all.search[!is.na(all.search$NCBI.ID),]
	if(nrow(all.search)==0){
		cat(paste("\t...No", clade.to.search, "sequences retrieved for any loci...\n",sep=" "))
	} else {
		cat(paste("\t...Retrieving rank level for", nrow(all.search), clade.to.search, "taxa...\n",sep=" "))
		if(nrow(all.search)<=500){
			tax.names<-NA
			attempt<-1
			while(is.na(tax.names) && attempt<=5){
				if(attempt>1){
					Sys.sleep(runif(1,min=1,max=10))
				}	
				try(tax.names<-ncbi_get_taxon_summary(id=all.search$NCBI.ID))
				attempt<-attempt+1
			}	
		} else {
			tax.names<-data.frame(matrix(ncol=3, nrow=nrow(all.search)))
			colnames(tax.names)<-c("uid","name","rank")
			for(p in 1:ceiling(nrow(all.search)/500)){
				if(p<ceiling(nrow(all.search)/500)){
					attempt<-1
					while(is.na(tax.names[((500*p)-499),]) && attempt<=5){
						if(attempt>1){
							Sys.sleep(runif(1,min=1,max=10))
						}	
						try(tax.names[((500*p)-499):(500*p),]<-ncbi_get_taxon_summary(id=all.search$NCBI.ID[((500*p)-499):(500*p)]))
						attempt<-attempt+1
					}	
				} else {
					attempt<-1
					while(is.na(tax.names[((500*p)-499),]) && attempt<=5){
						if(attempt>1){
							Sys.sleep(runif(1,min=1,max=10))
						}	
					try(tax.names[((500*p)-499):nrow(all.search),]<-ncbi_get_taxon_summary(id=all.search$NCBI.ID[((500*p)-499):nrow(all.search)]))
					attempt<-attempt+1
					}
				}
			}
		}
		non.sp<-tax.names$uid[tax.names$rank!="species"]
		tax.names<-tax.names[!tax.names$uid %in% non.sp,]
		all.search<-all.search[!all.search$NCBI.ID %in% non.sp,]
		all.search$Taxon<-tax.names$name[tax.names$uid %in% all.search$NCBI.ID]
		cat(paste("\t...Retrieving NCBI higher-level taxonomy for", nrow(all.search), clade.to.search, "taxa...\n",sep=" "))
		tax.output<-NULL
		attempt<-1
		if(nrow(all.search)<=500){
			attempt<-1
			while(is.null(tax.output) && attempt<=5){
				if(attempt>1){
					Sys.sleep(runif(1,min=1,max=500))
				}
				try(tax.output<-tax_name(all.search$Taxon,get=c("Species","Genus","Tribe","Subfamily","Family"),db="ncbi",verbose=FALSE))
				attempt<-attempt+1
			}		
		} else {
			for(p in 1:ceiling(nrow(all.search)/500)){
				if(p<ceiling(nrow(all.search)/500)){
					attempt<-1
					if(p==1){
						while(is.null(tax.output) && attempt<=5){
							if(attempt>1){
								Sys.sleep(runif(1,min=1,max=500))
							}
							try(tax.output<-tax_name(all.search$Taxon[((500*p)-499):(500*p)],get=c("Species","Genus","Tribe","Subfamily","Family"),db="ncbi",verbose=FALSE))
							attempt<-attempt+1
						}
					} else {
						attempt<-1
						while((nrow(tax.output)<(p*500)) && attempt<=5){
							if(attempt>1){
								Sys.sleep(runif(1,min=1,max=500))
							}	
							try(tax.output[((500*p)-499):(500*p),]<-tax_name(all.search$Taxon[((500*p)-499):(500*p)],get=c("Species","Genus","Tribe","Subfamily","Family"),db="ncbi",verbose=FALSE))
							attempt<-attempt+1
						}		
					}	
				} else {
					attempt<-1
					while((nrow(tax.output)<nrow(all.search)) && attempt<=5){
						if(attempt>1){
							Sys.sleep(runif(1,min=1,max=500))
						}
						try(tax.output[((500*p)-499):nrow(all.search),]<-tax_name(all.search$Taxon[((500*p)-499):nrow(all.search)],get=c("Species","Genus","Tribe","Subfamily","Family"),db="ncbi",verbose=FALSE))
						attempt<-attempt+1
					}
				}
			}
		}
		#return(tax.output)
		#the above line is useful if run into problems, such as when a taxon marked as "species" but lacks a proper taxon name due to it being an environmental sequence.
		#tax.output<-try(tax_name(all.search$Taxon,get=c("Species","Genus","Family","Order","Class","Subphylum","Phylum","Kingdom"),db="ncbi",verbose=FALSE))
		#this below should filter out things that have NA in species name in tax.output
		if(any(is.na(tax.output$Species))){
			tax.output<-tax.output[!is.na(tax.output$Species),]
			cat(paste("\t\t...Filtering accessions lacking proper names...\n",sep=" "))			
			cat(paste("\t...Retrieved NCBI data for", nrow(tax.output), clade.to.search, "taxa...\n",sep=" "))
		}
		if(all(tax.output$query==tax.output$Species)==TRUE){
			cat(paste("\t...All query names match species names...\n",sep=" "))
		}
		if(any(tax.output$query==tax.output$Species)==FALSE){
			cat(paste("\t...PROBLEM! All", clade.to.search, "query names DO NOT match species names...\n",sep=" "))
		}
		tax.output<-tax.output[,-c(1:2)]
		colnames(tax.output)[1]<-"Taxon"
		combined<-join(all.search,tax.output,by="Taxon",type="full")
		#This removes any taxa in all.search that lack Genus info
		combined<-combined[!is.na(combined$Genus),]
		combined<-arrange(combined,Family,Subfamily,Tribe,Genus,Taxon)
		return(combined)		
	}
}

write.acc.no.files<-function(df.w.info,path.for.files,clade,loci,append=FALSE){
	write.csv(df.w.info,file=paste(path.for.files,"/",clade,".combined.acc.nos.csv",sep=""),row.names=FALSE)
}

loci<-c("nuSSU","nuLSU","AbdA","Wg","LR","COI","EF1aF1","EF1aF2","ArgK","CAD","Top1","Ubx")
search.terms<-c("((18S[TITL] OR 16S[TITL] OR small[TITL]) AND ribosomal[TITL] NOT (internal[TITL] OR 58s[TITL] OR 5.8S[TITL] OR endonuclease[TITL] OR mitochondrial[TITL] OR mitochondrion[TITL] OR pseudogene [TITL] OR segmented set[PROP]","((28S[TITL] OR 26S[TITL] OR 25S[TITL] OR large[TITL]) AND ribosomal[TITL] NOT (internal[TITL] OR 58s[TITL] OR 5.8S[TITL] OR endonuclease[TITL] OR mitochondrial[TITL] OR mitochondrion[TITL] OR pseudogene [TITL] OR segmented set[PROP]","(AbdA[TITL] OR abd-A[TITL] OR abdominal-A[TITL] NOT (pseudogene [TITL] OR segmented set[PROP]","(Wg[TITL] wnt-1[TITL] OR wingless[TITL] NOT (pseudogene [TITL] OR segmented set[PROP]","(LR[TITL] OR LW Rh[TITL] OR (long[TITL] AND wavelength[TITL] AND rhodopsin[TITL]) OR (long[TITL] AND wave[TITL] AND opsin[TITL]) NOT (pseudogene [TITL] OR segmented set[PROP]","(COI[TITL] OR (cytochrome[TITL] AND oxidase[TITL] AND subunit[TITL] AND I[TITL]) OR (cytochrome[TITL] AND oxidase[TITL] AND subunit[TITL] AND 1[TITL]) NOT (pseudogene [TITL] OR segmented set[PROP]","(EF-1a-F1[TITL] OR EF1aF1[TITL] NOT (pseudogene [TITL] OR shotgun[TITL] OR scaffold[TITL] OR segmented set[PROP]","(EF-1a-F2[TITL] OR EF1aF2[TITL] NOT (pseudogene [TITL] OR shotgun[TITL] OR scaffold[TITL] OR segmented set[PROP]","(arginine kinase[TITL] ArgK[TITL] NOT (pseudogene [TITL] OR shotgun[TITL] OR scaffold[TITL] OR segmented set[PROP]","(carbomoylphosphate synthase[TITL] OR carbamoyl-phosphate synthetase[TITL] OR conserved ATPase domain protein isoform 1[TITL] OR CAD[TITL] NOT (pseudogene [TITL] OR shotgun[TITL] OR scaffold[TITL] OR segmented set[PROP]","(DNA topoisomerase 1[TITL] Top1[TITL] NOT (pseudogene [TITL] OR shotgun[TITL] OR scaffold[TITL] OR segmented set[PROP]","(ultrabithorax[TITL] Ubx[TITL] NOT (pseudogene [TITL] OR shotgun[TITL] OR scaffold[TITL] OR segmented set[PROP]")
excludes<-"OR AY233672 OR AY919019 OR AY233581 OR KC420256 OR KF908823 OR KF908825 OR DQ226038 OR EU439634 OR EU439628 OR FJ824312 OR KJ861361 OR DQ353252 OR EU561528 OR LN609174 OR LC025612 OR HQ440136 OR HQ440134 OR HQ440173 OR EU155441 OR EU155425 OR EU155429 OR EU155430 OR EU155431 OR EU155432 OR EU155434 OR EU155435 OR EU155436 OR EU155437 OR EU155438 OR EU155439 OR EU155442 OR DQ176297 OR NC_001566 OR DQ353402 OR DQ402089 OR EF012988 OR KP406340 OR KP406341 OR KP406342 OR KP406343 OR KP406344 OR KJ860035 OR EF013037 OR KJ859995 OR KJ859886 OR KF908822 OR AF147049 OR HQ853310 OR AY452149 OR AY452150 OR AY452148 OR AY265970 OR AY265967 OR AY265968 OR KJ855904 OR AY762648 OR JQ868225 OR JQ868081 OR EU143196 OR EU143194 OR EU143199 OR EU143192 OR EU143183 OR EU143184 OR EU143185 OR EU143203 OR EU143182 OR EU143206 OR EU143204 OR EU143193 OR EU143190 OR EU143186 OR EU143234 OR EU143191 OR EU143187 OR EU143188 OR EU143205 OR EU143197 OR EU143209 OR KJ443649 OR KJ443603 OR KJ443634 OR KJ443600 OR KJ443608 OR KJ443610 OR JX310584 OR JF759824 OR AF147043 OR AF147048 OR AF147051 OR FJ770201 OR FJ770206 OR FJ770203 OR FJ770204 OR HQ853327 OR DQ105549 OR AY519439 OR AY452158 OR AY519441 OR AY265964 OR AY265971 OR AF016015 OR AF016013 OR EU165430 OR AF016017 OR AY762645 OR AY762651 OR KJ855899 OR AB426047 OR AB668535 OR AB426043 OR HQ440138 OR HQ853328 OR JQ868069 OR JQ868076 OR EU143208 OR JQ617430 OR KJ443648 OR KJ443602 OR KJ443629 OR KJ443599 OR KJ443607 OR KJ443609 OR KF908830))"
clades<-"Formicidae"
path.for.files<-"SPECIFY_PATH"
sink(paste(path.for.files,"/",".acc.nos.txt",sep=""), append=TRUE, split=TRUE)

for(cl in 1:length(clades)){
	search.res<-get.accs.for.genes.in.clade(clade.to.search=clades[cl],locus.list=loci,search.terms=search.terms)
	if(!is.null(search.res)){
		write.acc.no.files(search.res,path.for.files,clades[cl],loci,append=FALSE)
	}
}
