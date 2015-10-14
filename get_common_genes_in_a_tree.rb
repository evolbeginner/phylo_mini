#! /bin/env ruby

require 'getoptlong'
require 'bio'


#####################################################################
tree_file = nil
gene_info_file = nil
genes_included_file = nil
target_taxon = nil
is_OTU_bootstrap = false

is_normal = false
is_rate = false
is_regular_norm = false
is_prop = false
is_ratio = false

target_node = nil
other_nodes = Array.new
gene_info = Hash.new
genes_included = Hash.new
all_genes = Hash.new

is_parent = false

ori_lengths = Array.new



#####################################################################
def get_genes_included(genes_included_file, genes_included=Hash.new)
  File.open(genes_included_file, 'r').each_line do |line|
    line.chomp!
    genes_included[line] = ""
  end
  return(genes_included)
end


def get_gene_info(gene_info_file, genes_included)
  gene_info = Hash.new
  File.open(gene_info_file, 'r').each_line do |line|
    line.chomp!
    line_array = line.split("\t")
    taxon_name = line_array[0]
    if ! genes_included.empty?
      line_array[1,line_array.size-1].each{|i| line_array.delete(i) if not genes_included.include?(i)}
    end
    gene_info[taxon_name] = line_array[1,line_array.size-1]
  end
  return(gene_info)
end


def get_descendents(tree, node, is_name)
  descendents = Array.new
  tree.descendents(node).each do |i|
    item = is_name ? i.name : i
    descendents.push(item)
  end
  return(descendents)
end


def get_parent_node(tree, node)
  adjacent_nodes = tree.adjacent_nodes(node)
  parent_node = adjacent_nodes[-1]
  return(parent_node)
end



def get_node_ancestral_genes(gene_info, tree, node, target_taxon, node_ancestral_genes, ancestral_genes={})
  node_ancestral_genes = []
  node_ancestral_genes_temp = Hash.new{|h,k|h[k]=[]}

  if node.name =~ /\w/
    if gene_info.include?(node.name)
      node_ancestral_genes = gene_info[target_taxon] & gene_info[node.name]
    else
      node_ancestral_genes = []
    end

  else
    tree.children(node).each do |sub_node|
      #next if sub_node.name == target_taxon
      if not ancestral_genes.include?(sub_node)
        ancestral_genes = get_node_ancestral_genes(gene_info, tree, sub_node, target_taxon, node_ancestral_genes, ancestral_genes)
      end

      descendents = get_descendents(tree, sub_node, true)
      if sub_node.name == target_taxon or descendents.include?(target_taxon)
        node_ancestral_genes_temp["target"].push(ancestral_genes[sub_node])
      else
        node_ancestral_genes_temp["non-target"].push(ancestral_genes[sub_node])
      end

      #p ancestral_genes
      node_ancestral_genes |= gene_info[target_taxon] & ancestral_genes[sub_node]
      #p node_ancestral_genes
      descendents = Hash.new
    end

    if node_ancestral_genes_temp.include?("target")
      node_ancestral_genes = node_ancestral_genes_temp["non-target"][0]
      node_descendents = get_descendents(tree, node, false)
      parent_nodes = tree.nodes - node_descendents
      parent_genes = Array.new
      parent_nodes.find_all{|i|i.name=~/\w/}.each do |i|
        next if not gene_info.include?(i.name)
        parent_genes |= gene_info[i.name]
      end
      node_ancestral_genes |= parent_genes
    else
      node_ancestral_genes = node_ancestral_genes_temp["non-target"].reduce([]){|sum, value| sum | value}
    end
  end

  ancestral_genes[node] = node_ancestral_genes
  return ancestral_genes
end



def get_minmax_bootstrap(tree, is_parent=false)
  bootstraps = Array.new
  tree.nodes.each do |node|
    if is_parent
      next if node.name =~ /\w/
    end
    parent_node = get_parent_node(tree, node)
    bootstrap = node.bootstrap - parent_node.bootstrap
    bootstraps << bootstrap
  end
  return(bootstraps.minmax)
end



#####################################################################
opts = GetoptLong.new(
  ['-t', '--tree', GetoptLong::REQUIRED_ARGUMENT],
  ['-g', '--gene', '--gene_table', '--gene_info', GetoptLong::REQUIRED_ARGUMENT],
  ['--genes_included', '--genes_included_file', GetoptLong::REQUIRED_ARGUMENT],
  ['--target', GetoptLong::REQUIRED_ARGUMENT],
  ['--OTU_bootstrap', GetoptLong::NO_ARGUMENT],
  ['--is_normal', GetoptLong::NO_ARGUMENT],
  ['--is_rate', GetoptLong::NO_ARGUMENT],
  ['--is_regular_norm', GetoptLong::NO_ARGUMENT],
  ['--is_prop', GetoptLong::NO_ARGUMENT],
  ['--ratio', '--is_ratio', GetoptLong::NO_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '-t', '--tree'
      tree_file = value
    when '-g', '--gene', '--gene_table', '--gene_info'
      gene_info_file = value
    when '--genes_included', '--genes_included_file'
      genes_included_file = value
    when '--target'
      target_taxon = value
    when '--OTU_bootstrap'
      is_OTU_bootstrap = true
    when '--is_normal'
      is_normal = true
    when '--is_rate'
      is_rate = true
    when '--is_regular_norm'
      is_regular_norm = true
    when '--is_prop'
      is_prop = true
    when '--ratio', '--is_ratio'
      is_ratio = true
  end
end


[tree_file, gene_info_file].each do |file|
  if file.nil?
    raise "File #{file} has to be given! Exiting ......"
  end
end


#####################################################################
genes_included = get_genes_included(genes_included_file, genes_included) if ! genes_included_file.nil?
gene_info = get_gene_info(gene_info_file, genes_included)


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


ancestral_genes = Hash.new{|h,k|h[k]=[]}


other_nodes.each do |node|
  if is_OTU_bootstrap
    if node.name =~ /\w/
      #p tree.adjacent_nodes(node)
      node.bootstrap = "#"+(gene_info[target_taxon] & gene_info[node.name]).size.to_s
    end
  end
  node_ancestral_genes = Array.new
  ancestral_genes = get_node_ancestral_genes(gene_info, tree, node, target_taxon, node_ancestral_genes, ancestral_genes)
  node.bootstrap = ancestral_genes[node].size
end


################################################################################
is_parent=true if is_prop

min_bootstrap, max_bootstrap = get_minmax_bootstrap(tree, is_parent)
tree.nodes.each do |node|
  parent_node = get_parent_node(tree, node)
  if is_normal
    length = node.bootstrap - parent_node.bootstrap
  elsif is_rate
    if tree.get_edge(node, parent_node).distance.to_f != 0.0
      length = (node.bootstrap - parent_node.bootstrap)/tree.get_edge(node, parent_node).distance.to_f
    end
  elsif is_regular_norm
    length = (node.bootstrap - parent_node.bootstrap - min_bootstrap)/(max_bootstrap - min_bootstrap).to_f
  elsif is_prop
    length = (node.bootstrap - parent_node.bootstrap)/[node.bootstrap, parent_node.bootstrap].max.to_f
  end
  ori_lengths << length
  tree.get_edge(node, parent_node).distance = length
end


ori_lengths << nil
ori_lengths.compact!
ori_length_min = ori_lengths.compact.min
tree.nodes.each do |node|
  parent_node = get_parent_node(tree, node)
  ori_length = tree.get_edge(node, parent_node).distance
  if ori_length
    tree.get_edge(node, parent_node).distance = ori_length + 1.2 * (ori_lengths.min).abs
  end
end



################################################################################
gene_info.values.each do |i|
  i.each do |j|
    all_genes[j] = ""
  end
end

tree.nodes.each do |node|
  node.bootstrap = (node.bootstrap.to_f/all_genes.size*100).round(1).to_s+"%" if is_ratio
  if node.name =~ /\w/
    node.bootstrap = ""
  end
end


tree.options[:bootstrap_style]=:traditional
output = tree.output_newick
output.gsub!(/[\n\s]/, '')

puts output


