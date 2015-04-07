#! /usr/bin/perl
use strict;
use Getopt::Long;


#########################################################################
my ($newick_file, $bootstrap_cutoff, $branch_length_cutoff);
my ($is_paml_label);
$bootstrap_cutoff = 0;
$branch_length_cutoff=0;

GetOptions(
	'in|newick_file=s'	=>	\$newick_file,
	'bootstrap=s'		=>	\$bootstrap_cutoff,
	'length|branch_length=s'=>	\$branch_length_cutoff,
	'label|paml_label!'	=>	\$is_paml_label,
) || die "illegal param!\n";


#########################################################################
open (my $IN, '<', "$newick_file");
while(<$IN>){
	next if $_ !~ /\w/;
	$_ =~ s/(?<=\))([0-9]+(\.\d+)?)/{$1 >= $bootstrap_cutoff ? $1 : ''}/eg;
	$_ =~ s/\:([^:;, ()]+)/{$1 >= $branch_length_cutoff ? ':'.$1 : ''}/eg;
	if ($is_paml_label){
	    my $label_counter;
	    $label_counter=0;
	    $_ =~ s/(?=[,)])/{"#".$label_counter++}/eg if $is_paml_label;
	}
	print;
}


