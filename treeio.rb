#! /usr/bin/env ruby


#########################################################
# A script to manipulate the trees in Newick format using BioRuby
# Last updated: 2018-04-09
# Author: Sishuo Wang from Haiwei Luo Lab at Chinese University of Hong Kong (CUHK)
# E-mail: sishuowang@hotmail.ca
# Notes: 1. This script is actively being updated. So its usage may be frequently changed. Please be aware of the version. Sorry for the inconvenience.
#        2. The script together with its associated scripts is distributed under CC-BY 0.0, in the hope that it will be useful, but WITHOUT ANY WARRANTY.
#        3. Please cite the following paper when using it. Goto, N. et al. BioRuby: bioinformatics software for the Ruby programming language. Bioinformatics 26, 2617–2619 (2010).


#########################################################
dir = File.dirname($0)
$: << File.join(dir, 'treeRubyLib')


#########################################################
require 'getoptlong'
require 'rubystats'

require 'bio'
require 'tree'


#########################################################
tree_file = nil
is_output_otu = false
bootstrap_min = -10000
bootstrap_style = :traditional
name_bootstrap_style = nil
bl_times = 1
bl_scale = nil
bs_times = 1
bl_min = 0
is_bl_one = false
bootstrap_select = Array.new
is_bs_rev = false
is_digital = false
is_output_parent_child = false
is_dist = false
is_remove_poly = false
is_print_node_name = false


#########################################################
def output_output_parent_child_rela(tree)
  nodeCount = 0
  tree.nodes.each do |node|
    if node.name !~ /\w/
      nodeCount += 1
      node.name = (nodeCount).to_s
    end
  end
  tree.each_node do |node|
    parent_name = tree.parent(node).nil? ? nil : tree.parent(node).name
    puts [node.name, parent_name, tree.children(node).map{|child| child.name}].join("\t")
  end
end


def nameBootstrap(tree, name_bootstrap_style)
  if not name_bootstrap_style.nil?
    tree.internal_nodes.each_with_index do |internal_node, index|
      taxa = tree.tips(internal_node).map{|tip|tip.name.gsub(' ','_')}.join('-')
      case name_bootstrap_style
        when :taxa
          internal_node.bootstrap = taxa
        when :numeric
          internal_node.bootstrap = '!' + (index+1).to_s
      end
    end
  end
  return(tree)
end


def splitBootstrap(tree, bootstrap_select)
  tree.internal_nodes.each do |node|
    if node.name =~ /[\/]/
      if bootstrap_select
        bs = Array.new
        bootstrap_select.each do |i|
          bs << node.name.split('/')[i-1]
        end
        node.bootstrap = bs.join('/')
        node.name = nil
      end
    end
  end
  return(tree)
end


def revBootstrap(tree)
  tree.internal_nodes.each do |node|
    node.bootstrap = -(node.bootstrap.to_f)
  end
  return(tree)
end


def output_otu(tree)
  #tree.each_node do |node|
  #  puts node.bootstrap
  #end; exit
  tree.allTips.each do |tip|
    puts tip.name
  end
end


#########################################################
opts = GetoptLong.new(
  ['-i', '-t', GetoptLong::REQUIRED_ARGUMENT],
  ['-b', GetoptLong::REQUIRED_ARGUMENT],
  ['--output_otu', '--output_OTU', GetoptLong::NO_ARGUMENT],
  ['--no_b', GetoptLong::NO_ARGUMENT],
  ['--name_bootstrap', GetoptLong::REQUIRED_ARGUMENT],
  ['--bl_times', GetoptLong::REQUIRED_ARGUMENT],
  ['--bl_scale', GetoptLong::REQUIRED_ARGUMENT],
  ['--bs_times', GetoptLong::REQUIRED_ARGUMENT],
  ['--bl_min', GetoptLong::REQUIRED_ARGUMENT],
  ['--bl_one', GetoptLong::NO_ARGUMENT],
  ['--bootstrap_select', '--bs_select', GetoptLong::REQUIRED_ARGUMENT],
  ['--bs_rev', '--rev_bs', GetoptLong::NO_ARGUMENT],
  ['--digit', '--digital', GetoptLong::NO_ARGUMENT],
  ['--output_parent_child', GetoptLong::NO_ARGUMENT],
  ['--dist', GetoptLong::NO_ARGUMENT],
  ['--remove_poly', GetoptLong::NO_ARGUMENT],
  ['--print_node_name', GetoptLong::NO_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^-(i|t)$/
      tree_file = value
    when /^--output_(otu|OTU)$/
      is_output_otu = true
    when /^-b$/
      bootstrap_min = value.to_f
    when /^--no_b$/
      bootstrap_style = :disabled
    when /^--name_bootstrap$/
      name_bootstrap_style = value.to_sym
    when /--bl_times/
      bl_times = value.to_f
    when /--bl_scale/
      bl_scale = value
    when /--bs_times/
      bs_times = value.to_f
    when /--bl_min/
      bl_min = value.to_f
    when /--bl_one/
      is_bl_one = true
    when /^--(bootstrap_select|bs_select)$/
      bootstrap_select = value.split(',').map{|i|i.to_i}
    when /^--(bs_rev|rev_bs)$/
      is_bs_rev = true
    when '--digit', '--digital'
      is_digital = true
    when '--output_parent_child'
      is_output_parent_child = true
    when '--dist'
      is_dist = true
    when '--remove_poly'
      is_remove_poly = true
    when '--print_node_name'
      is_print_node_name = true
  end
end


#########################################################
#treeio = Bio::FlatFile.open(Bio::Newick, tree_file)
trees = getTreeObjs(tree_file)


#########################################################
def get_norm_distr(bl_scale)
  mu, sigma = bl_scale.split(',').map{|i|i.to_f}
  norm_distr = Rubystats::NormalDistribution.new(mu, sigma)
  return(norm_distr)
end


# randomly generated beta-distributed params to scale the branch length
norm_distr = get_norm_distr(bl_scale) if not bl_scale.nil? 


#########################################################
trees.each do |tree|
  tree.options[:bootstrap_style] = bootstrap_style
  if is_output_otu
    output_otu(tree)
    exit
  end

  output_output_parent_child_rela(tree) if is_output_parent_child

  tree = nameBootstrap(tree, name_bootstrap_style)

  tree = splitBootstrap(tree, bootstrap_select)

  tree = revBootstrap(tree) if is_bs_rev


  #########################################################
  tree.each_node do |node|
    #node.bootstrap_string = '>10<15' unless node.isTip?(tree)
    #puts tree.distance(tree.root, node) if is_dist
    if is_print_node_name
      puts tree.twoTaxaNodeName(node).map{|i|i.gsub(' ', '_')}.join('|') unless node.isTip?(tree); next
    end

    if is_dist
      next if tree.root == node
      puts tree.distance(tree.parent(node), node)
    end

    if node.respond_to?(:bootstrap) and not node.bootstrap.nil?
      #if not node.bootstrap.is_a?(String)
      if not node.bootstrap.is_a?(Numeric)
        node.bootstrap = (node.bootstrap * bs_times).to_i
        if node.bootstrap < bootstrap_min
          node.bootstrap = nil
        end
      end
    end

    if is_remove_poly
      if tree.children(node).size >= 3
        poly_node_names = tree.tips(node).map{|i|i.name.gsub(' ', '_')}
        p poly_node_names
        puts poly_node_names.shuffle[0,poly_node_names.size-2].join("\n")
      end
    end
  end


  tree.each_edge do |node0, node1, edge|
    next if edge.distance.to_s !~ /\d/
    edge.distance = edge.distance + bl_min if edge.distance < bl_min
    edge.distance = edge.distance * bl_times
    edge.distance = edge.distance * [norm_distr.rng(), 0.1].max if not bl_scale.nil? # scale bl with a norm distr
    edge.distance = 1 if is_bl_one
    if is_digital
      edge.distance = "%0.9f" % edge.distance
    end
    if bl_times == 0
      edge.distance = nil
    end
    #p edge.distance
  end


  output = tree.output_newick
  output.gsub!(/[\n\s]/, '')

  puts output
end


