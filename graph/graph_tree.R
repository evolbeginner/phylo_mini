#! /bin/env Rscript


#######################################################
# plot a phylogenetic tree with "ggtree" in the format of PDF
# rerooting supported
# written by Sishuo Wang (sishuowang@hotmail.ca)


#######################################################
#suppressWarnings(suppressMessages(library("ggtree", quietly = T)))
library('ggtree')
library('getopt')
library('ape')
library('phangorn')


#######################################################
root = c()
xlim_add = c(0, 1)


#######################################################
spec <- matrix(
	c("infile",  "i", 1, "character", "the infile",
    	"outfile", "o", 1, "character", "the outfile",
    	"no_bl", "", 0, "logical", "no branch length",
    	"root", "r", 1, "character", "the root(s)",
    	"midpoint", "m", 0, "logical", "midpoint root the tree",
    	"xlim", "x", 1, "character", "xlim"
	),
	byrow=TRUE, ncol=5
)

opt <- getopt(spec=spec)

infile = opt$infile
outfile = opt$outfile


#######################################################
if (infile == '-'){
	#tree<-read.tree(text="(Organism1.006G249400.1:10.03977,(Organism2.022118m:0.01337,(Organism3.J34265.1:0.00284,Organism4.G02633.1:0.00468)0.51:0.0104):0.02469);")}
	;
}else{
	tree<-read.tree(infile)
}

if ("root" %in% names(opt) ){
	roots = unlist(strsplit(opt$root, ","))
	tree <- root(tree, outgroup=roots, edgelabel = TRUE)
}

if ("no_bl" %in% names(opt)){
	tree <- ggtree(tree, branch.length="none")
}

if ("midpoint" %in% names(opt)){
	getRoot(tree)
	tree <- midpoint(tree)
}

if ("xlim" %in% names(opt)){
	xlim_add = opt$xlim
}


#######################################################
pdf(outfile)

#print(xlim_add)

#ggtree(tree) + geom_text2(aes(subset=!isTip, label=label), size=2, hjust = -0.1) + geom_text2(aes(subset=isTip, label=label), hjust = -0.3) + xlim(xlim_add)
ggtree(tree) + geom_text2(aes(subset=!isTip, label=label), size=1, hjust = 0) + geom_text2(aes(subset=isTip, label=label), hjust = 10) + xlim(xlim_add)

dev.off()


