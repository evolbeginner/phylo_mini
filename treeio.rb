#! /bin/env ruby

require 'bio'

#########################################################
tree_file = ARGV[0]

treeio = Bio::FlatFile.open(Bio::Newick, tree_file)
tree = treeio.next_entry.tree
tree.options[:bootstrap_style] = :traditional

output = tree.output_newick
output.gsub!(/[\n\s]/, '')

puts output


