#! /bin/env perl

use strict;
use 5.010;
use Getopt::Long;
use read_gene_seq_from_paml;


#######################################################################
my (%seq_info, %kaks_info, $treefile);

my $create_codemlctl = "/home/sswang/tools/self_bao_cun/PAML/create_codemlctl.pl";

my ($input, @targets, $seqtype, $is_help);
#my ($create_codemlctl, $input, $treefile, $outfile, $ncatG, $seqtype);


#######################################################################
GetOptions(
    'i|input=s'     =>	\$input,    # paml format
    't|target=s'    =>  \@targets,
    'seqtype=s'     =>  \$seqtype,
	'help|h!'		=>	\$is_help,
) || die "illegal params! Exiting ......";

&usage() if $is_help;

my $seq_info_href = read_gene_seq_from_paml::read_paml($input);
%seq_info = %$seq_info_href;

$treefile = $input . ".tree";
&create_triplet_tree_for_paml([keys(%seq_info)], $treefile);

&run_codeml();

my $kaks_info_href = &get_kaks_treeview("mlc", $seqtype);
%kaks_info = %$kaks_info_href;

print join("\t", map {my $target=$_; $kaks_info{$target} if exists $kaks_info{$target}} @targets) . "\n";


#######################################################################
sub run_codeml{
	no strict;
	my ($output_codemlctl_dir, $NSsites);
	$output_codemlctl_dir="./";
	$NSsites = 0;
	$model = 2;
	system"\"$create_codemlctl\" -seqfile=\"$input\" -treefile=\"$treefile\" -outfile=\"$outfile\" -model=\"$model\" -NSsites=\"$NSsites\" -ncatG=\"$ncatG\" -seqtype=\"$seqtype\" -output_codemlctl_dir=\"$output_codemlctl_dir\"";
	`codeml`;
}


sub get_kaks_treeview{
    my %kaks_info;
    my ($mlc, $seqtype) = @_;
    my $sign;
    if ($seqtype == "2"){
        $sign = "tree length";
    }
    else{
        $sign = 'w ratios as labels for TreeView:';
    }
    open (my $IN, '<', $mlc) || die "mlc cannot be opened!";
    while(my $line = <$IN>){
        #print "!!\n" if $line =~ /$sign/;
        chomp($line);
        if ($line =~ /$sign/){
            if ($seqtype == "2"){
                <$IN>; <$IN>;
            }
            while($line = <$IN>){
                chomp($line);
                next if $line =~ /^$/;
                #(AT1G28150 #0.4208 , Bra032879_LF #0.2316 , Bra030075_MF2 #0.3077 );
                #(Bo7g116340_MF1: 0.20085, Bo1g007900_LF: 0.12526, AT4G33720: 0.10842);
                my @a;
                my @line = split(/,/, $line);
                if ($seqtype == "2"){
                    map {$_=~s/[()#; ]//g} @line;
                }
                else{
                    map {$_=~s/[()#;]//g} @line;
                }
                if ($seqtype == "2"){
                    @a = map {(split(":"))} @line;
                }
                else{
                    @a = map {(split)} @line;
                }
                foreach (0..$#a-1){
                    $kaks_info{$a[$_]} = $a[$_+1];
                }
                last;
            }
        }
    }
    return (\%kaks_info);
}


sub get_kaks_pairwise{
    my (@taxa, %kaks_info);
    my ($mlc) = @_;
    open (my $IN, '<', $mlc) || die "mlc cannot be opened!";
    my $sign = 'Nei & Gojobori 1986.';
    #my $sign = 'Nei & Gojobori 1986. dN/dS (dN, dS)';
    while(my $line = <$IN>){
        chomp($line);
        if ($line =~ /^$sign/){
            my ($is_start, $taxa_counter);
            $taxa_counter = 0;
            while(my $line = <$IN>){
                chomp($line);
                if ($line =~ /^$/){
                    $is_start = 1;
                    next;
                }

                if ($is_start){
                    my (%kaks_info_each, $counter);
                    $counter = 1;
                    $taxa_counter += 1;
                    my @line = split(/\s+/, $line);
                    my $taxon = $line[0];
                    push @taxa, $taxon;
                    foreach my $ele (1..$#line){
                        $line[$ele] =~ s/[()]//g;
                        given($ele){
                            when ($ele%3 == 1){
                                $kaks_info{$taxa_counter}{$counter}{'kaks'} = $line[$ele];
                            }
                            when ($ele%3 == 2){
                                $kaks_info{$taxa_counter}{$counter}{'ka'} = $line[$ele];
                            }
                            when ($ele%3 == 0){
                                $kaks_info{$taxa_counter}{$counter}{'ks'} = $line[$ele];
                            }
                        }
                        if ($ele%3 == 0){
                            #print join("\t", $counter, $taxa_counter, $kaks_info{$taxa_counter}{$counter}{kaks})."\n";
                            $counter += 1;
                        }                        
                    }
                }

            }
        }
    }

    foreach my $key1 (keys %kaks_info){
        foreach my $key2 (keys %{$kaks_info{$key1}}){
            print $key1, $key2, "\n";
            print $kaks_info{$key1}{$key2}{'kaks'}."\n";
        }
    }
    close ($IN);
}


sub create_triplet_tree_for_paml{
    my ($seq_titles_aref, $treefile) = @_;
    open(my $OUT, '>', $treefile) || die "treefile $treefile cannot be created!";
    print $OUT "(";
    my @new_taxa = map {$seq_titles_aref->[$_]."\#$_"} 0..scalar(@$seq_titles_aref)-1;
    print $OUT join(",", @new_taxa);
    print $OUT ");\n";
    close $OUT;
}


sub usage{
    print "The usage of $0 is:";
    exit;
}


