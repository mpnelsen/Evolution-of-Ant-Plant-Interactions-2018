#This python file can be used to BLAST sequences against a local database
#The file trim.hybrid.family2.r is called from within and used to parse output

import re
import os
import csv
import time
import shutil
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from cogent import LoadSeqs, DNA
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, Reference

date="06july2016"
loci = ["nuSSU", "nuLSU", "AbdA", "COI", "LR", "Wg", "EF1aF1", "EF1aF2", "ArgK", "CAD", "Top1", "Ubx"]
folderpath="PATH
megablast="/PATHTO/blast-2.2.26/bin/megablast"
blastdb="/PATHTO/antblastdb/combined_refseqs_parsed_UPPER.db"
trim_r="/PATHTO/trim.hybrid.family2.r"
clade="Formicidae"

for locus in loci:
	#perform local blast to check if seqs actually from organisms of interest
	subprocess.call(args="{0} -i {1}{2}/{3}/{3}_Renamed_Parsed.fasta -o {1}{2}/{3}/{3}_blast.txt -W 8 -r 2 -q -3 -G 5 -E 2 -v 1 -b 1 -m 8 -d {4}".format(megablast,folderpath,clade,locus,blastdb), shell=True);
	shutil.copy("{0}{1}/{2}/{2}_blast.txt".format(folderpath,clade,locus),"{0}{1}/{2}/{2}_blast_for_edit.txt".format(folderpath,clade,locus));

	#use some perl to edit file
	subprocess.call(args="perl -i -pe 's/^[#].*\n.*/\1/g' {0}{1}/{2}/{2}_blast_for_edit.txt".format(folderpath,clade,locus), shell=True);
	subprocess.call(args="perl -i -pe 's/^(\S+\t\S+).*/$1/g' {0}{1}/{2}/{2}_blast_for_edit.txt".format(folderpath,clade,locus), shell=True);
	#don't need the \001 line here if running perl natively (ie, not from w/in python)
	subprocess.call(args="perl -pi -e 's/\001//g' {0}{1}/{2}/{2}_blast_for_edit.txt".format(folderpath,clade,locus), shell=True);
	subprocess.call(args="uniq {0}{1}/{2}/{2}_blast_for_edit.txt > {0}{1}/{2}/{2}_blast_for_edit_clean.txt".format(folderpath,clade,locus), shell=True);


subprocess.call(args="Rscript {0}".format(trim_r),shell=True);

#this then pulls out those that blast correctly...
for locus in loci:
	aln = LoadSeqs('{0}{1}/{2}/{2}_Renamed_Parsed.fasta'.format(folderpath,clade,locus), moltype=DNA, aligned=False)
	keepers_first_round=open("{0}{1}/{2}/{2}.taxa.blastok.txt".format(folderpath,clade,locus))
	have=keepers_first_round.read().split()
	have_align=aln.takeSeqs(have)
	have_align.writeToFile('{0}{1}/{2}/{2}_blastok.fasta'.format(folderpath,clade,locus));
