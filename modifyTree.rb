#! /bin/env ruby


require 'getoptlong'
require 'bio'


#############################################################
tree_file = nil


#############################################################
def get_parent_node(tree, node)
  adjacent_nodes = tree.adjacent_nodes(node)
  parent_node = adjacent_nodes[-1]
  return(parent_node)
end


#############################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      tree_file = value
  end
end


#############################################################
treeio = Bio::FlatFile.open(Bio::Newick, tree_file)
tree = treeio.next_entry.tree


tree.nodes.each do |node|
  parent_node = get_parent_node(tree, node)
  next if tree.get_edge(node, parent_node).distance.nil?
  next if node.name !~ /\w/
  node.name.sub!(/^[^|]+\|/, '')
  #tree.get_edge(node, parent_node).distance *= 10
end


output = tree.output_newick
output.gsub!(/[\n\s]/, '')

puts output

