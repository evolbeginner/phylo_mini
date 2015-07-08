#! /bin/env ruby

require 'getoptlong'
require 'bio'


#####################################################################
tree_file = nil
gene_info_file = nil
target_taxon = nil
is_OTU_bootstrap = false

target_node = nil
other_nodes = Array.new
gene_info = Hash.new


#####################################################################
def get_gene_info(gene_info_file)
  gene_info = Hash.new
  File.open(gene_info_file, 'r').each_line do |line|
    line.chomp!
    line_array = line.split("\t")
    taxon_name = line_array[0]
    gene_info[taxon_name] = line_array[1,line_array.size-1]
  end
  return(gene_info)
end


#####################################################################
opts = GetoptLong.new(
  ['-t', '--tree', GetoptLong::REQUIRED_ARGUMENT],
  ['-g', '--gene', '--gene_table', '--gene_info', GetoptLong::REQUIRED_ARGUMENT],
  ['--target', GetoptLong::REQUIRED_ARGUMENT],
  ['--OTU_bootstrap', GetoptLong::NO_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '-t', '--tree'
      tree_file = value
    when '-g', '--gene', '--gene_table', '--gene_info'
      gene_info_file = value
    when '--target'
      target_taxon = value
    when '--OTU_bootstrap'
      is_OTU_bootstrap = true
  end
end


[tree_file, gene_info_file].each do |file|
  if file.nil?
    raise "File #{file} has to be given! Exiting ......"
  end
end


#####################################################################
gene_info = get_gene_info(gene_info_file)


treeio = Bio::FlatFile.open(Bio::Newick, tree_file)
tree = treeio.next_entry.tree
tree.options[:bootstrap_style] = :disabled

tree.each_node {|x| other_nodes.push(x)}
nodes = tree.nodes

nodes.each do |node|
  if node.name == target_taxon
    target_node = node
  end
end

other_nodes.each do |other_node|
  next if other_node.name !~ /\w/
  #puts other_node
  #p tree.lowest_common_ancestor(target_node, other_node)
  #puts (gene_info[target_taxon] & gene_info[other_node.name]).join("\t")
  #puts "#######"
end


#puts other_nodes.map{|i|i.name}.join("\t")


other_nodes.each do |node|
  #tree.descendents(i).each{|x| print x.name.to_s}; puts
  if is_OTU_bootstrap
    if node.name =~ /\w/
      p tree.adjacent_nodes(node)
      node.bootstrap = "#"+(gene_info[target_taxon] & gene_info[node.name]).size.to_s
    end
  end
  next if node.name =~ /\w/
  node_ancestral_genes = Array.new
  #puts ["node:", node.name].join("\t")
  tree.descendents(node).each do |sub_node|
    next if sub_node.name !~ /\w/
    next if sub_node.name == target_taxon
    node_ancestral_genes |= gene_info[target_taxon] & gene_info[sub_node.name]
    #p tree.children(node)
  end
  #puts tree.descendents(node).find_all{|i|i.name=~/\w/}.map{|i|i.name}.join("\t")
  #puts node_ancestral_genes.join("\t")
  node.bootstrap = node_ancestral_genes.size
  #print tree.adjacent_nodes(i); print "\t"
end


output = tree.output_newick
output.gsub!(/[\n\s]/, '')

puts output

#sub_tree.each_node{|x| puts x.name}
#puts tree.distance(tree.nodes[0,1].to_a)
#puts tree.parent(a[0],a[2])


