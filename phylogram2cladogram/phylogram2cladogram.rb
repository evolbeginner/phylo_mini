#! /usr/bin/env ruby


#########################################################
DIR = File.dirname(__FILE__)
$: << File.join(DIR, 'lib')


#########################################################
require 'getoptlong'

require 'tree'


#########################################################
infile = nil

addDists = Hash.new


#########################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
)


#########################################################
opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
  end
end


#########################################################
tree = getTreeObjs(infile)[0]
tree2 = getTreeObjs(infile)[0]

root = tree.root
tree.each_edge do |node0, node1, edge|
  edge.distance = 1
end
tree2.each_edge do |node0, node1, edge|
  edge.distance = 1
end


tree.each_node do |node|
  next if node == root
  dist = tree.distance(node,root)
  sister = tree.sister(node)
  tips_dists = tree.tips(node).map{|i|tree.distance(i,root)}
  sister_tips_dists = tree.tips(sister).map{|i|tree.distance(i,root)}
  parent = tree.parent(node)
  if tips_dists.max < sister_tips_dists.max
    taxaName = tree.twoTaxaNodeName(node)
    addDists[taxaName] = sister_tips_dists.max - tips_dists.max
    #p [node, tree.sister(node), dist, tips_dists.max, edge.distance]
  end
end


#########################################################
tree = tree2
tree.each_node do |node|
  taxaName = tree.twoTaxaNodeName(node)
  parent = tree.parent(node)
  if addDists.include?(taxaName)
    edge = tree.get_edge(node, parent)
    edge.distance += addDists[taxaName]
  end
end


#########################################################
puts tree.cleanNewick()


