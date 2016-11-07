#! /bin/env Rscript

library("ape")

args <- commandArgs(trailingOnly = TRUE)

infile = args[1]
type = args[2]

intree = read.tree(args[1])

plot(intree, type=type, edge.width=5, cex=2)

