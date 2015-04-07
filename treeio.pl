#! /usr/bin/perl
use Bio::TreeIO;
use Bio::Tree::Tree;


$treefile = $ARGV[0];
$outfile  = $ARGV[1];

die "Must provide treefile as the 1st param." unless defined $treefile;

my $TreeStream  = new Bio::TreeIO(
					'-format' => do {defined $ARGV[2] ? $ARGV[2] : 'nexus'},
					'-file'   => $treefile
					);

my $tree_obj = $TreeStream->next_tree();
my $outTree = Bio::TreeIO->new(
	-file => ">$outfile",
	-format => do {defined $ARGV[3] ? $ARGV[3] : 'newick'},
#	-internal_node_id => 'bootstrap'
);
	$outTree->write_tree($tree_obj);
