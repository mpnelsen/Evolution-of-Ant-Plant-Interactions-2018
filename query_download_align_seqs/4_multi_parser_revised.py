#python file to parse sequences from downloaded NCBI files.

#! /usr/bin/env python
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os, sys
loci=["nuSSU","nuLSU","AbdA","Wg","LR","COI","EF1aF1","EF1aF2","ArgK","CAD","Top1","Ubx"]
#pathtofolders="/Users/matthewnelsen/Desktop"
pathtofolders="PATH"
locus_genenames={}
locus_genenames["nuSSU"]=["18S ribosomal RNA","small subunit ribosomal RNA","18S small subunit ribosomal RNA", "small subunit 18S ribosomal RNA", "sequence contains 18S rRNA gene", "18S rRNA", "nuclear encoded 18S ribosomal RNA, small subunit", "nuclear encoded small subunit ribosomal RNA", "small subunit (18S) ribosomal RNA", "18S ribosomal RNA", "16S small subunit ribosomal RNA", "18S small ribosomal RNA subunit", "nuclear small subunit ribosomal RNA", "contains 18S ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2, and 28S ribosomal RNA", "16S ribosomal RNA", "SSU ribosomal RNA","contains 18S ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2, and 26S ribosomal RNA", "contains 18S ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2, and 28S ribosomal RNA", "small subunit ribosomal RNA","18S SSU ribosomal RNA", "18S small subunit ribosomal RNA", "ribosomal RNA small subunit","contains 18S ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2, and 28S ribosomal RNA","small ribosomal RNA subunit RNA","SSU ribosomal RNA","18S rRNA","nuclear 18S ribosomal RNA","18S ribosomal RNA gene","put. 18S ribosomal RNA","18S ribosomal RNA, small subunit"]
locus_genenames["nuLSU"]=["28S ribosomal RNA","large subunit ribosomal RNA","28S large subunit ribosomal RNA", "large subunit 28S ribosomal RNA", "sequence contains 28S rRNA gene", "28S rRNA", "nuclear encoded 28S ribosomal RNA, large subunit", "nuclear encoded large subunit ribosomal RNA", "large subunit (28S) ribosomal RNA", "28S ribosomal RNA", "25S large subunit ribosomal RNA", "28S large ribosomal RNA subunit", "nuclear large subunit ribosomal RNA", "25S ribosomal RNA", "LSU ribosomal RNA", "26S rRNA", "26S ribosomal RNA", "26S large subunit ribosomal RNA", "large subninit ribosomal RNA", "large subunit ribosomal RNA", "contains 18S ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2, and 28S ribosomal RNA", "contains 18S ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2, and 26S ribosomal RNA", "sequence contains 28S rRNA gene", "28S large subunit ribosomal","large subunit ribosomal RNA","sequence contains 28S rRNA gene","contains internal transcribed spacer 2 and 28S ribosomal RNA","sequence contains 18S rRNA gene, ITS1, 5.8S rRNA gene, ITS2, 28S rRNA gene", "28S small subunit ribosomal RNA","put. 26S ribosomal RNA","5.8S ribosomal RNA and large subunit ribosomal RNA"]

def parse(file):
	for record in SeqIO.parse(open(file,"r"),"genbank"):
		accession=record.id
		species=record.annotations['organism']
		species2=species.replace(" ","_")
		species3=species2.replace(".","")
		if locus=="nuSSU" or locus=="nuLSU":
			print locus
			for i, feature in enumerate (record.features):
				if locus=="nuSSU" or locus=="nuLSU":
				#print('{0}_genenames'.format(locus))
				#genenames=["18S ribosomal RNA","small subunit ribosomal RNA","18S small subunit ribosomal RNA", "small subunit 18S ribosomal RNA", "sequence contains 18S rRNA gene", "18S rRNA", "nuclear encoded 18S ribosomal RNA, small subunit", "nuclear encoded small subunit ribosomal RNA", "small subunit (18S) ribosomal RNA", "18S ribosomal RNA", "16S small subunit ribosomal RNA", "18S small ribosomal RNA subunit", "nuclear small subunit ribosomal RNA"]
					if any(feature.type=='rRNA' for i, feature in enumerate (record.features)):
						if feature.type=='rRNA':
							if 'product' in feature.qualifiers:
								geneofinterest=feature.qualifiers['product'][0]
								seq = feature.extract(record.seq)
								if geneofinterest in locus_genenames[locus]:
									print '>%s %s\n%s' % (species3, geneofinterest, seq)
									parsedseqfile.write('>{0}\n{1}\n'.format(accession,seq))
									acclistfile.write('%s\t%s\n' % (species3,accession))
									break
							else:
								if 'note' in feature.qualifiers:
									geneofinterest=feature.qualifiers['note'][0]
									seq = feature.extract(record.seq)
									if geneofinterest in locus_genenames[locus]:
										print '>%s %s\n%s' % (species3, geneofinterest, seq)
										parsedseqfile.write('>{0}\n{1}\n'.format(accession,seq))
										acclistfile.write('%s\t%s\n' % (species3, accession))
										break		
					else:
						if any(feature.type=='misc_feature' for i, feature in enumerate (record.features)):
							if feature.type=='misc_feature':
								if 'note' in feature.qualifiers:
									geneofinterest=feature.qualifiers['note'][0]
									seq = feature.extract(record.seq)
									if geneofinterest in locus_genenames[locus]:
										print '>%s %s\n%s' % (species3, geneofinterest, seq)
										parsedseqfile.write('>{0}\n{1}\n'.format(accession,seq))
										acclistfile.write('%s\t%s\n' % (species3, accession))
										break		
						else:
							if any(feature.type=='misc_RNA' for i, feature in enumerate (record.features)):
								if feature.type=='misc_RNA':
									if 'product' in feature.qualifiers:
										geneofinterest=feature.qualifiers['product'][0]
										seq = feature.extract(record.seq)
										if geneofinterest in locus_genenames[locus]:
											print '>%s %s\n%s' % (species3, geneofinterest, seq)
											parsedseqfile.write('>{0}\n{1}\n'.format(accession,seq))
											acclistfile.write('%s\t%s\n' % (species3, accession))
											break		
									else:
										if 'note' in feature.qualifiers:
											geneofinterest=feature.qualifiers['note'][0]
											seq = feature.extract(record.seq)
											if geneofinterest in locus_genenames[locus]:
												print '>%s %s\n%s' % (species3, geneofinterest, seq)
												parsedseqfile.write('>{0}\n{1}\n'.format(accession,seq))
												acclistfile.write('%s\t%s\n' % (species3,accession))
												break
							else:
								if any(feature.type=='precursor_RNA' for i, feature in enumerate (record.features)):
									if feature.type=='precursor_RNA':
										if 'product' in feature.qualifiers:
											geneofinterest=feature.qualifiers['product'][0]
											seq = feature.extract(record.seq)
											if geneofinterest in locus_genenames[locus]:
												print '>%s %s\n%s' % (species3, geneofinterest, seq)
												parsedseqfile.write('>{0}\n{1}\n'.format(accession,seq))
												acclistfile.write('%s\t%s\n' % (species3,accession))
												break
										else:
											if 'note' in feature.qualifiers:
												geneofinterest=feature.qualifiers['note'][0]
												seq = feature.extract(record.seq)
												if geneofinterest in locus_genenames[locus]:
													print '>%s %s\n%s' % (species3, geneofinterest, seq)
													parsedseqfile.write('>{0}\n{1}\n'.format(accession,seq))
													acclistfile.write('%s\t%s\n' % (species3, accession))
													break		
		if locus=="AbdA" or locus=="Wg" or locus=="LR" or locus=="COI" or locus=="EF1aF1" or locus=="EF1aF2" or locus=="ArgK" or locus=="CAD" or locus=="Top1" or locus=="Ubx":
			if any(feature.type=='CDS' for i, feature in enumerate (record.features)):
				for i, feature in enumerate (record.features):
					if feature.type=='CDS':
						if 'product' in feature.qualifiers:
							seq = feature.extract(record.seq)
							print '>%s %s\n%s' % (species3, locus, seq)
							parsedcdsseqfile.write('>{0}\n{1}\n'.format(accession,seq))
							cdsacclistfile.write('%s\t%s\n' % (species3, accession))
							if 'translation' in feature.qualifiers:
								protein=feature.qualifiers['translation'][0]
								print '>%s %s\n%s' % (species3, locus, protein)
								parsedprotfile.write('>{0}\n{1}\n'.format(accession,protein))
								protacclistfile.write('%s\t%s\n' % (species3, accession))
								break
			else:
				if any(feature.type=='gene' for i, feature in enumerate (record.features)):
					for i, feature in enumerate (record.features):
						if feature.type=='gene':
							if 'gene' in feature.qualifiers:
								seq = feature.extract(record.seq)
								print '>%s %s\n%s' % (species3, locus, seq)
								parsedgeneseqfile.write('>{0}\n{1}\n'.format(accession,seq))
								geneacclistfile.write('%s\t%s\n' % (species3, accession))
								break
				else:
					if any(feature.type=='misc_feature' for i, feature in enumerate (record.features)):
						for i, feature in enumerate (record.features):
							if feature.type=='misc_feature':
								seq = feature.extract(record.seq)
								print '>%s %s\n%s' % (species3, locus, seq)
								parsedmiscseqfile.write('>{0}\n{1}\n'.format(accession,seq))
								miscacclistfile.write('%s\t%s\n' % (species3, accession))
								break		

for locus in loci:
	if locus=="AbdA" or locus=="Wg" or locus=="LR" or locus=="COI" or locus=="EF1aF1" or locus=="EF1aF2" or locus=="ArgK" or locus=="CAD" or locus=="Top1" or locus=="Ubx":
		parsedcdsseqfile = open("{0}/{1}/{1}_CDS_Parsed.fasta".format(pathtofolders,locus),"wa")
		cdsacclistfile = open("{0}/{1}/{1}_CDS_accession_list.txt".format(pathtofolders,locus),"wa")
		parsedgeneseqfile = open("{0}/{1}/{1}_GENE_Parsed.fasta".format(pathtofolders,locus),"wa")
		geneacclistfile = open("{0}/{1}/{1}_GENE_accession_list.txt".format(pathtofolders,locus),"wa")
		parsedmiscseqfile = open("{0}/{1}/{1}_MISC_Parsed.fasta".format(pathtofolders,locus),"wa")
		miscacclistfile = open("{0}/{1}/{1}_MISC_accession_list.txt".format(pathtofolders,locus),"wa")
		parsedprotfile = open("{0}/{1}/{1}_AA_Parsed.fasta".format(pathtofolders,locus),"wa")
		protacclistfile = open("{0}/{1}/{1}_AA_accession_list.txt".format(pathtofolders,locus),"wa")
		parse("{0}/{1}/{1}.gb".format(pathtofolders,locus))
		parsedcdsseqfile.close()
		cdsacclistfile.close()
		parsedgeneseqfile.close()
		geneacclistfile.close()
		parsedmiscseqfile.close()
		miscacclistfile.close()
		parsedprotfile.close()
		protacclistfile.close()
	if locus=="nuSSU" or locus=="nuLSU":
		parsedseqfile = open("{0}/{1}/{1}_Parsed.fasta".format(pathtofolders,locus),"wa")
		acclistfile = open("{0}/{1}/{1}_accession_list.txt".format(pathtofolders,locus),"wa")
		parse("{0}/{1}/{1}.gb".format(pathtofolders,locus))
		parsedseqfile.close()
		acclistfile.close();
		
