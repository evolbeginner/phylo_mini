#! /bin/env Rscript

library("phangorn")
#library("phytools")

args <- commandArgs(trailingOnly = TRUE)

infile = args[1]
outfile = args[2]

intree = read.tree(args[1])
outtree = midpoint(intree)
#outtree = midpoint.root(intree)

write.tree(outtree, file=outfile)

