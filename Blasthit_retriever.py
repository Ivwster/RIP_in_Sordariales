#!/usr/bin/env python3

import os
import sys
import pysam
import re
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Filters TE blastn hits based on coverage and percent identity and fetches fasta files of TE blastn hits")

parser.add_argument("-b","--blast",help="The blast output file (outfmt 6)",type=str)
parser.add_argument("-s","--strain",help="the name of the strain",type=str)
parser.add_argument("-t","--TE",help="TE library Fasta file",type=str)
parser.add_argument("-g","--genome",help="TE library Fasta file",type=str)
parser.add_argument("-c","--cov",help="coverage threshold. Default = 0",type=int,default=0)
parser.add_argument("-l","--len",help="Length threshold. Default = 0",type=int,default=0)
parser.add_argument("-p","--pident",help="percent identity threshold. Default = 0",type=int,default=0)
#parser.add_argument("-m","--mode",help="'all' repeats or 'split' based on superfamily, default = all",default="all",type=str)
parser.add_argument("-o","--output",help="Output directory, where the output files are created", type=str)

args = parser.parse_args()

if len(sys.argv) <= 1:
	parser.print_help()
	sys.exit()

if args.TE is None:
    parser.print_help()
    sys.exit("Missing TE fasta input file")

if args.strain is None:
    parser.print_help()
    sys.exit("Missing strain input")

if args.blast is None:
    parser.print_help()
    sys.exit("Missing blast table file")

if args.genome is None:
    parser.print_help()
    sys.exit("Missing genome fasta file")

if args.output is None:
    parser.print_help()
    sys.exit("Missing output directory")



TE_file=args.TE
strain=args.strain
blast_file=args.blast
genome=args.genome
cov=args.cov
perc_ident=float(args.pident)
out=args.output
len_thresh=args.len


# Retrieves fasta file information from specific chromosome and positions
def retriever(fasta,chrom, start, end):
	seq = pysam.FastaFile(fasta)
	cur_seq = seq.fetch(chrom, start-1, end-1)
    
	return cur_seq


blast_tab=[]
library_families=[]

with open (blast_file) as blast:
     
     for hit in blast:
          
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = hit.strip().split("\t")

        blast_tab.append((qseqid, sseqid, float(pident), int(length), mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore))
    
with open(TE_file) as library:
    
    for line in library:

        #print(line) 
        if ">" in line:
         
            family,clas = line.strip().split("#")

            fam = family.replace(family[0],"",1) 

            library_families.append((fam,clas))


with open (f"{out}{strain}_blastn_copies_filtered_{cov}-cov_{perc_ident}-pident_{len_thresh}-len.out","w") as out_1:
    with open (f"{out}{strain}_blastn_copies_filtered_{cov}-cov_{perc_ident}-pident_{len_thresh}-len.fasta","w") as out_2:

        for hit in blast_tab:
     
            for family in library_families:
        
        
                if hit[0].lower() == f"{family[0].lower()}#{family[1].lower()}":
            
                    query_seq = pysam.FastaFile(TE_file)

                    cur_seq = query_seq.fetch(f"{family[0]}#{family[1]}")

                    #print(hit[0],"\t",hit[1],"\t",hit[3]/len(cur_seq),"\t",cov/100,"\t",hit[2],"\t",perc_ident)
                    if hit[3]/len(cur_seq) >= cov/100 and hit[3] >= len_thresh and hit[2] >= perc_ident:
          
                
                    
                        out_1.write(f"{hit[0]}\t{hit[1]}\t{hit[2]}\t{hit[3]}\t{hit[4]}\t{hit[5]}\t{hit[6]}\t{hit[7]}\t{hit[8]}\t{hit[9]}\t{hit[10]}\t{hit[11]}\n")
        
                
             
                        hit_seq = retriever(genome,hit[1],min(int(hit[8]),int(hit[9])),max(int(hit[8]),int(hit[9])))

                        out_2.write(f">{hit[0]}_{hit[1]}_{hit[8]}-{hit[9]}\n{hit_seq}\n")