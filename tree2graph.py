#! /bin/env python


import sys
import getopt

from Bio import Phylo
import pylab


#####################################################################
def get_label(leaf):
    return leaf.name


#####################################################################
dpi = 400
out_fmt = "pdf"


#####################################################################
try:
    opts, args = getopt.getopt(
        sys.argv[1:],
        "i:o:",
        ["in=", "out=", "outfmt=", "out_fmt=", "dpi="],
    )
except getopt.GetoptError:
    print "Illegal params!"
    sys.exit()


for op, value in opts:
    if op == "-i" or op == "--in":
        infile = value
    elif op == "-o" or op == "--out":
        outfile = value
    elif op == "--outfmt" or op == "--out_fmt":
        out_fmt = value
    elif op == "--dpi":
        dpi = int(value)


#####################################################################
tree = Phylo.read(infile, 'newick')
tree.ladderize()
Phylo.draw(tree, label_func=get_label, do_show=False)
pylab.axis('off')
pylab.savefig(outfile, format=out_fmt, bbox_inches='tight', dpi=dpi)


