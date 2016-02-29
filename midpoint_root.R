#! /bin/env Rscript


library("phytools")
args <- commandArgs(trailingOnly = TRUE)

infile = args[1]
outfile = args[2]

intree = read.newick(args[1])
outtree = midpoint.root(intree)

write.tree(outtree, file=outfile)

