#! /bin/env ruby

require 'getoptlong'
require 'bio'


#####################################################################
tree_file = nil
gene_info_file = nil
genes_included_file = nil
genes_excluded_file = nil
target_taxon = nil
cluster_file = nil
cluster_type = "orthomcl"
is_OTU_bootstrap = false
length_divided_by = nil
length_Jian_by = nil
length_ChengFang_by = nil
outgroup_taxa = Array.new

is_normal = false
is_rate = false
is_regular_norm = false
is_prop = false
is_ratio = false
is_no_change_length = false
is_2ndMin = false

target_node = nil
other_nodes = Array.new
gene_info = Hash.new
genes_included = Hash.new
genes_excluded = Hash.new
all_genes = Hash.new
clusters = Hash.new
gene_cluster_info = Hash.new

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


def get_gene_info(gene_info_file, genes_included, genes_excluded)
  gene_info = Hash.new
  File.open(gene_info_file, 'r').each_line do |line|
    line.chomp!
    line_array = line.split("\t")
    taxon_name = line_array[0]
    if ! genes_included.empty?
      line_array[1,line_array.size-1].each{|i| line_array.delete(i) if not genes_included.include?(i)}
    end
    if ! genes_excluded.empty?
      line_array[1,line_array.size-1].each{|i| line_array.delete(i) if genes_excluded.include?(i)}
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


def read_cluster_file(cluster_file, type="orthomcl")
  in_fh = File.open(cluster_file, 'r')
  clusters = Hash.new{|h,k|h[k]=[]}
  gene_cluster_info = Hash.new
  count = 0
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    count += 1
    if type == "orthomcl"
      line_arr.each do |i|
        i=$1 if i =~ /[^|]+\|(.+)/
        clusters[count] << i
        gene_cluster_info[i] = count
      end
    elsif type == "silix"
      line_arr.each do |i|
        clusters[line_arr[0]] << line_arr[1]
        gene_cluster_info[line_arr[1]] = line_arr[0]
      end
    end
  end
  in_fh.close
  return([clusters, gene_cluster_info])
end


def obtain_full_gene_cluster_info(gene_cluster_info, clusters, gene_info, target_taxon)
  count = clusters.size - 1
  gene_info[target_taxon].each do |gene|
    if ! gene_cluster_info.include?(gene)
      count += 1
      gene_cluster_info[gene] = count
      #p gene_cluster_info.size
    end
  end
  return(gene_cluster_info)
end


#####################################################################
opts = GetoptLong.new(
  ['-t', '--tree', GetoptLong::REQUIRED_ARGUMENT],
  ['-g', '--gene', '--gene_table', '--gene_info', GetoptLong::REQUIRED_ARGUMENT],
  ['--genes_included', '--genes_included_file', GetoptLong::REQUIRED_ARGUMENT],
  ['--genes_excluded', '--genes_excluded_file', GetoptLong::REQUIRED_ARGUMENT],
  ['--target', GetoptLong::REQUIRED_ARGUMENT],
  ['--cluster_file', '--clstr', '--cluster', '--clstr_file', GetoptLong::REQUIRED_ARGUMENT],
  ['--cluster_type', '--clstr_type', GetoptLong::REQUIRED_ARGUMENT],
  ['--OTU_bootstrap', GetoptLong::NO_ARGUMENT],
  ['--is_normal', GetoptLong::NO_ARGUMENT],
  ['--is_rate', GetoptLong::NO_ARGUMENT],
  ['--is_regular_norm', GetoptLong::NO_ARGUMENT],
  ['--is_prop', GetoptLong::NO_ARGUMENT],
  ['--is_no_change_length', '--no_change_length', GetoptLong::NO_ARGUMENT],
  ['--ratio', '--is_ratio', GetoptLong::NO_ARGUMENT],
  ['--length_divided_by', GetoptLong::REQUIRED_ARGUMENT],
  ['--length_Jian_by', GetoptLong::REQUIRED_ARGUMENT],
  ['--length_ChengFang_by', GetoptLong::REQUIRED_ARGUMENT],
  ['--outgroup', GetoptLong::REQUIRED_ARGUMENT],
  ['--is_2ndMin', '--2ndMin', GetoptLong::NO_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '-t', '--tree'
      tree_file = value
    when '-g', '--gene', '--gene_table', '--gene_info'
      gene_info_file = value
    when '--genes_included', '--genes_included_file'
      genes_included_file = value
    when '--genes_excluded', '--genes_excluded_file'
      genes_excluded_file = value
    when '--target'
      target_taxon = value
    when '--cluster_file', '--clstr', '--cluster', 'clstr_file'
      cluster_file = value
    when '--cluster_type', '--clstr_type'
      cluster_type = value
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
    when '--is_no_change_length', '--no_change_length'
      is_no_change_length = true
    when '--ratio', '--is_ratio'
      is_ratio = true
    when '--length_divided_by'
      length_divided_by = value.to_f
    when '--length_Jian_by'
      length_Jian_by = value.to_f
    when '--length_ChengFang_by'
      length_ChengFang_by = value.to_f
    when '--outgroup'
      outgroup_taxa << value
    when '--2ndMin'
      is_2ndMin = true
  end
end


[tree_file, gene_info_file].each do |file|
  if file.nil?
    raise "File #{file} has to be given! Exiting ......"
  end
end


if not cluster_file.nil?
  if cluster_file =~ /silix$/
    cluster_type = "silix"
  end
end


#####################################################################
genes_included = get_genes_included(genes_included_file, genes_included) if ! genes_included_file.nil?
genes_excluded = get_genes_excluded(genes_excluded_file, genes_excluded) if ! genes_excluded_file.nil?

gene_info = get_gene_info(gene_info_file, genes_included, genes_excluded)

clusters, gene_cluster_info = read_cluster_file(cluster_file, cluster_type) if ! cluster_file.nil?
gene_cluster_info = obtain_full_gene_cluster_info(gene_cluster_info, clusters, gene_info, target_taxon)


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
  if ! clusters.empty?
    families = Hash.new
    ancestral_genes[node].map{|i|families[gene_cluster_info[i]] = ""}
    node.bootstrap = families.size
  else
    node.bootstrap = ancestral_genes[node].size
  end
end


################################################################################
is_parent=true if is_prop

min_bootstrap, max_bootstrap = get_minmax_bootstrap(tree, is_parent)
tree.nodes.each do |node|
  parent_node = get_parent_node(tree, node)
  next if outgroup_taxa.include?(node.name) 
  if is_normal
    length = node.bootstrap - parent_node.bootstrap
  elsif is_regular_norm
    length = (node.bootstrap - parent_node.bootstrap - min_bootstrap)/(max_bootstrap - min_bootstrap).to_f
  elsif is_prop
    length = (node.bootstrap - parent_node.bootstrap)/[node.bootstrap, parent_node.bootstrap].max.to_f
  elsif is_no_change_length
    next if tree.get_edge(node, parent_node).distance.nil?
    #tree.get_edge(node, parent_node).distance *= 10
    next
  end
  if is_rate
    if tree.get_edge(node, parent_node).distance.to_f != 0.0
      length /= tree.get_edge(node, parent_node).distance.to_f
    end
  end
  ori_lengths << length
  tree.get_edge(node, parent_node).distance = length
end


ori_lengths << nil
ori_lengths.compact!
if is_2ndMin
  ori_length_min = (ori_lengths.sort)[1]
else
  ori_length_min = (ori_lengths.sort)[0]
end

tree.nodes.each do |node|
  parent_node = get_parent_node(tree, node)
  ori_length = tree.get_edge(node, parent_node).distance
  if ori_length and not is_no_change_length
    tree.get_edge(node, parent_node).distance = ori_length + 1.2 * (ori_length_min).abs
    #puts tree.get_edge(node, parent_node).distance
    tree.get_edge(node, parent_node).distance /= length_divided_by if ! length_divided_by.nil?
    tree.get_edge(node, parent_node).distance -= length_Jian_by if ! length_Jian_by.nil?
    tree.get_edge(node, parent_node).distance **= length_ChengFang_by if ! length_ChengFang_by.nil?
  end
end



################################################################################
gene_info.values.each do |i|
  if clusters.empty?
    i.each do |j|
      all_genes[j] = ""
    end
  else
    i.each do |j|
      all_genes[gene_cluster_info[j]] = ""
    end
  end
end


tree.options[:bootstrap_style]=:traditional

tree.nodes.each do |node|
  if node.name =~ /\w/
    prop = nil
    if is_ratio
      if ! clusters.empty?
        prop = (node.bootstrap.to_f/all_genes.size*100).round(1).to_s+"%" if is_ratio
      else
        prop = (gene_info[node.name].size.to_f/all_genes.size*100).round(1).to_s+"%" if is_ratio
      end
    end
    node.name = prop.to_s + "-" + node.name
  end
end


tree.nodes.each do |node|
  node.bootstrap = (node.bootstrap.to_f/all_genes.size*100).round(1).to_s+"%" if is_ratio
  if node.name =~ /\w/
    node.bootstrap = ""
  end
end

output = tree.output_newick
output.gsub!(/[\n\s]/, '')


puts output



