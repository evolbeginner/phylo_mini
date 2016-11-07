#! /bin/env perl

use Bio::TreeIO;
use strict;
use 5.010;
use Getopt::Long;


#############################################################################
my ($tree_file, $outfile);


#############################################################################
GetOptions(
        'tree=s'    =>  \$tree_file,
        'out=s' =>  \$outfile,
) || die "parameters error!";


#############################################################################
my $in = new Bio::TreeIO(-file => $tree_file,
                         -format => 'newick');
my $out = new Bio::TreeIO(-file => ">$outfile",
                          -format => 'svggraph',
                          -font_size => "50");
while( my $tree = $in->next_tree ) {
    $out->write_tree($tree);
}


