#python file to fetch sequences

import sys
import time
import math
import os
from Bio import Entrez
from Bio import SeqIO


Entrez.email = "email"

loci = ["nuSSU", "nuLSU", "AbdA", "COI", "LR", "Wg", "EF1aF1", "EF1aF2", "ArgK", "CAD", "Top1", "Ubx"]
folderpath="PATHFORFILES"

print loci
batchsize=200
for locus in loci:
    os.mkdir("{0}/{1}".format(folderpath,locus))
    accnos = []
    try:
        for line in open('{0}/{1}.acc.nos.txt'.format(folderpath,locus), 'r'):
            accnos.append(line.strip('\n'))
            #os.mkdir("{0}/{1}".format(folderpath,locus))
            #os.mkdir("{0}/{1}/Raw_Seqs".format(folderpath,locus))		
        print "accnos: ", accnos
        #problemfile = open("{0}/{1}/{1}.seq.retrieval.problems.txt".format(folderpath, locus), 'w')
        #accnosstr=",".join(map(str,accnos))
        #from https://www.biostars.org/p/66921/
        #from https://www.biostars.org/p/63506/
        #get gi's for accession numbers (useful for epost, as that does not take accessions)
        iters=int(math.ceil(float(len(accnos))/batchsize))
        output = open("/{0}/{1}/{1}.gb".format(folderpath,locus), 'wa')
        for it in range(0,iters):
            gi_hand = Entrez.esearch(db="nucleotide",term= " ".join(accnos[it*batchsize:((it+1)*batchsize)]), retmax=batchsize)
            gis=Entrez.read(gi_hand)['IdList']
			#sometimes it has a hard time getting everything and i get a service temporarily unavailable note - by printing here  and sleeping below, it seems to slow it down enough
            print gis
            gi_res_hand= Entrez.epost(db="nucleotide", id=",".join(gis))
            print gi_res_hand
            gi_results = Entrez.read(gi_res_hand)
            print gi_results
            webEnv=gi_results["WebEnv"]
            queryKey=gi_results["QueryKey"]
            for gi in range(0,len(gis),batchsize):
                giHand=Entrez.efetch(db="nucleotide",rettype="gb",retmode="text",retstart=gi,retmax=batchsize,webenv=webEnv,query_key=queryKey)
                output.write(giHand.read())
                time.sleep(2)
    except:
        print "{0} NOT FOUND! CONTINUING TO NEXT LOCUS...".format(locus)

output.close();
