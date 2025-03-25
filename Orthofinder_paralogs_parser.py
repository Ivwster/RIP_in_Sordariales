#!/usr/bin/env python3

import os
import sys
import pysam
import re
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Duplications per species parser")

parser.add_argument("-g","--orthologues",help="The ortologues tsv file",type=str)
parser.add_argument("-d","--duplications",help="The duplications tsv file",type=str)
parser.add_argument("-s","--strain",help="the name of the text file that lists what strains to include",type=str)
parser.add_argument("-f","--fasta",help="Genome Fasta path",type=str)
parser.add_argument("-a","--annotation",help="Genome annotation path",type=str)
parser.add_argument("-o","--output",help="Output directory, where the output files are created", type=str)

args = parser.parse_args()

if len(sys.argv) <= 1:
	parser.print_help()
	sys.exit()

if args.orthologues is None:
    parser.print_help()
    sys.exit("Missing Orthologues.tsv file")

if args.duplications is None:
    parser.print_help()
    sys.exit("Missing Duplications.tsv file")

if args.strain is None:
    parser.print_help()
    sys.exit("Missing strain input file")

if args.fasta is None:
    parser.print_help()
    sys.exit("Missing genome fasta path")

if args.annotation is None:
    parser.print_help()
    sys.exit("Missing genome gene annotation input")

if args.output is None:
    parser.print_help()
    sys.exit("Missing output directory")


ortho_file=args.orthologues
dups_file=args.duplications
strains_file=args.strain
fasta_path=args.fasta
annotation_path=args.annotation
output=args.output


# Retrieves fasta file information from specific chromosome and positions
def retriever(fasta,chrom, start, end):
	seq = pysam.FastaFile(fasta)
	cur_seq = seq.fetch(chrom, start-1, end-1)
    
	return cur_seq


## Make path if it does not already exist
os.makedirs(output, exist_ok=True)

workingdir = f"{output}"
## make workingdir
os.makedirs(workingdir,exist_ok=True)

with open (strains_file) as strains:
    
    for cur_strain in strains:
        cur_strain=cur_strain.strip()

        if cur_strain == "neuint-8807":
            cur_strain_id = "8807pb"
        elif cur_strain == "neumet-10397":
            cur_strain_id = "10397pb"
        elif cur_strain == "neusit-5941":
            cur_strain_id = "5941pb"
        elif cur_strain == "neutre-9045":
            cur_strain_id = "9045pb"
        elif cur_strain == "podans-137":
            cur_strain_id= "PaWa137m"
        elif cur_strain == "podcom-139":
            cur_strain_id = "PcWa139m"
        elif cur_strain == "podbel-112042":
            cur_strain_id = "QC761"
        elif cur_strain == "podpse-415":
            cur_strain_id = "QC762"
        elif cur_strain == "sormac-2":
            cur_strain_id = "SMAC4"
        else:
            #cur_strain_id = cur_strain
            cur_strain_id = re.sub("\\-(.*)","",cur_strain.strip(" "))


        annotation_list=[]
        with open (f"{annotation_path}{cur_strain}.gff") as cur_annotation:
            
            for line in cur_annotation:
                if "#" in line: 
                    pass
                else:
                    lines = line.strip().split('\t')
                    annotation_list.append((lines))
                    

        ortho_dict = {}

        with open (ortho_file) as ortho:
            

            for i, line in enumerate(ortho):
                if i > 0:
                    
                    lis = line.strip().split('\t')
                    
                    for gene in lis:
                        
                        

                        if cur_strain_id.lower() in gene.lower():
                            

                            if len(gene.split(',')) > 1:
                                
                                
                                if cur_strain == "neuint-8807" or cur_strain == "neumet-10397" or cur_strain == "neusit-5941" or cur_strain == "neutre-9045":
                                    
                                    genes = gene.split(',')

                                    for j, g in enumerate(genes):
                                        genes[j] = re.sub("-mRNA-1","",re.sub(".*(?=gene)","",genes[j]))


                                    ortho_dict[lis[0]] = genes

                                elif cur_strain == "podpse-415" or cur_strain == "podcom-139" or cur_strain == "podbel-112042" or cur_strain == "podans-137" or cur_strain == "sormac-2":

                                    genes = gene.split(',')

                                    for j, g in enumerate(genes):
                                        genes[j] = re.sub(" ","",genes[j])
                                    
                                    ortho_dict[lis[0]] = genes
                                    
                                else:
                                    genes = gene.split(',')
                                    genes = genes
                                    for j, g in enumerate(genes):

                                        if genes[j][0].isupper() == True or genes[j][1].isupper() == True:
                                            genes[j] = genes[j].lower()
                                            genes[j] = re.sub(f"{cur_strain_id}_","",genes[j])
                                            
                                        else:
                                            genes[j] = re.sub("_t1","",re.sub(f"{cur_strain}_","",genes[j]))
                                    
                                    ortho_dict[lis[0]] = genes
        unique_genes=[]                        
        with open (f"{workingdir}{cur_strain}_paralogs_sequences.fasta", "w") as out:

            for d in ortho_dict:
                for gene in ortho_dict[d]:
                    gene = gene.strip(" ")
                    for an in annotation_list:
                        
                        if an[2] == "gene":
                            
                            
                            if f"{gene};" in an[8]:
                    
                                cur_seq = retriever(f"{fasta_path}{cur_strain}.faa",an[0],int(an[3]),int(an[4]))

                                if f"{cur_strain}_{d}_{gene}_{an[0]}_{an[3]}-{an[4]}" in unique_genes:
                                    pass
                                else:
                                    unique_genes.append(f"{cur_strain}_{d}_{gene}_{an[0]}_{an[3]}-{an[4]}")
                                    out.write(f">{cur_strain}_{d}_{gene}_{an[0]}_{an[3]}-{an[4]}\n{cur_seq}\n")
                            
                            elif f"ID={gene}" == an[8]:

                                cur_seq = retriever(f"{fasta_path}{cur_strain}.faa",an[0],int(an[3]),int(an[4]))
                                
                                if f"{cur_strain}_{d}_{gene}_{an[0]}_{an[3]}-{an[4]}" in unique_genes:
                                    pass

                                else:
                                    unique_genes.append(f"{cur_strain}_{d}_{gene}_{an[0]}_{an[3]}-{an[4]}")
                                    out.write(f">{cur_strain}_{d}_{gene}_{an[0]}_{an[3]}-{an[4]}\n{cur_seq}\n")
                                
        print(f"done with {cur_strain}")






