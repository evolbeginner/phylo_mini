#! /bin/env Rscript

library(ape)

args <- commandArgs(trailingOnly = TRUE)

t=read.tree(args[1])

t_unrooted = unroot(t)

cat(write.tree(t_unrooted))
