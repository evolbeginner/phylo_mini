library(ape)

args <- commandArgs(trailingOnly = TRUE)
infile = args[1]

x = read.dna(infile,format="fasta")
#dist.dna(x, model="RAW")
dist.dna(x, model="K80")
