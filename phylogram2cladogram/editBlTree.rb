#! /usr/bin/env ruby


###############################################################
DIR = File.dirname(__FILE__)
$: << File.join(DIR, 'lib')


###############################################################
require 'getoptlong'

require 'tree'


###############################################################
infile = nil
subtree_files = Array.new
zoom = 1

taxaNameIncluded = Array.new


###############################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--subtree', GetoptLong::REQUIRED_ARGUMENT],
  ['--zoom', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '--subtree'
      subtree_files << value.split(',')
    when '--zoom'
      zoom = value.to_f
  end
end


subtree_files.flatten!


###############################################################
subtree_files.each do |subtree_file|
  trees = getTreeObjs(subtree_file)
  trees.each do |tree|
    root = tree.root
    tree.each_node do |node|
      taxaNodeName = tree.twoTaxaNodeName(node)
      if node == root and taxaNodeName.size == 1
        taxaNodeName = tree.twoTaxaNodeName(tree.children(root)[0])
      end
      taxaNameIncluded << taxaNodeName
    end
  end
end


###############################################################
tree = getTreeObjs(infile)[0]
tree.each_node do |node|
  taxaNodeName = tree.twoTaxaNodeName(node)
  if taxaNameIncluded.include?(taxaNodeName)
    #taxaNodes = tree.twoTaxaNode(node)
    edge = tree.get_edge(node, tree.parent(node))
    edge.distance *= zoom
  end
end


###############################################################
puts tree.cleanNewick()


