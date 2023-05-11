#! /usr/bin/env ruby


####################################################
$DIR = File.dirname($0)
$: << File.join($DIR, 'lib')


####################################################
require 'getoptlong'
require 'bio'


require 'tree'


####################################################
infile = nil
tree_file = nil


####################################################
def getRoot(infile)
  taxa = Array.new
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    taxa << line
  end
  in_fh.close
  taxa.uniq.sort!
  return(taxa)
end


def getNewRootedTrees(trees, root_taxa)
  taxa2node = Hash.new
  trees.each do |tree|
    tree.nodes.each do |node|
      #next if node.isTip?(species_tree)
      taxa = tree.tips(node).map{|i|i.name}.sort
      taxa2node[taxa] = node
    end
    new_root_taxa = taxa2node.keys.select{|taxa|(taxa&root_taxa)==root_taxa}.sort_by{|taxa|(taxa-root_taxa).size}
    raise "Root taxa not found in the tree" if new_root_taxa.empty?
    new_root = taxa2node[new_root_taxa]
    tree.root = new_root
  end
  return(trees)
end


def outputTrees(trees)
  trees.each do |tree|
    puts tree.cleanNewick()
  end
end


####################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['-t', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when /^-i$/
      infile = value
    when /^-t$/
      tree_file = value
  end
end


####################################################
trees = getTreeObjs(tree_file, 10)

root_taxa = getRoot(infile)

trees = getNewRootedTrees(trees, root_taxa)

outputTrees(trees)


