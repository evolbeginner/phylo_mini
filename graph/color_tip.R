#! /bin/env Rscript


########################################################
options(warn=-1)
suppressMessages(library(getopt))
suppressMessages(library(ggtree))
suppressMessages(library(phytools))


########################################################
library(getopt)
library(ggtree)
library(phytools)


########################################################
font_size = 10
layout = 'rectangular'
is_midpoint = F
is_branch_length = 'branch.length'


########################################################
command=matrix(c( 
    'help', 'h', 0, 'logical',
    'tree', 't', 2, 'character',
    'annot', 'a', 2, 'character',
    'output', 'o', 2, 'character',
    'layout', 'l', 2, 'character',
    'midpoint', 'm', 0, 'logical',
    'no_bl', '', 0, 'logical',
    'size', 's', 2, 'double'),
    byrow=T, ncol=4
)
args=getopt(command)


########################################################
treefile <- args$tree
annot <- args$annot
output <- args$output

if(is.null(args[["layout"]])){
	;
}else{
	layout <- args$layout
}

if(! is.null(args[["size"]])){
	font_size <- args$size
}

if(! is.null(args[["midpoint"]])){
	is_midpoint <- T
}

if(! is.null(args[["no_bl"]])){
	is_branch_length <- "none"
}


########################################################
tree <- read.tree(treefile)
if(is_midpoint){
	tree <- midpoint.root(tree)
}
df <- read.table(annot, header=T)

p <- ggtree(tree, layout=layout, branch.length=is_branch_length, size=font_size/5)

#p <- p %<+% df + geom_point(aes(shape=isTip, color=isTip), size=0.5)
p <- p %<+% df + geom_tiplab(aes(color=cat), size=font_size) 
p <- p + theme(legend.position='none')

pdf(output)
p
dev.off()


