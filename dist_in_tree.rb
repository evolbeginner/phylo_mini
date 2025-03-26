#! /usr/bin/env ruby


###########################################
# to find distance btwn pairs of OTUs in a tree


###########################################
require 'getoptlong'
require 'bio'
require 'bio-nwk'

require 'util'


###########################################
treefile = nil
infile = nil

pairs = Array.new
distances = Array.new


###########################################
opts = GetoptLong.new(
  ['-t', GetoptLong::REQUIRED_ARGUMENT],
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '-t'
      treefile = value
    when '-i'
      infile = value
  end
end


###########################################
tree = getTreeObjs(treefile)[0]

in_fh = File.open(infile, 'r')
in_fh.each do |line|
  line.chomp!
  pairs << line.split("\t")
end
in_fh.close

pairs.map!{|a|a.map!{|i|i.gsub('_', ' ')}}

pairs.each do |genes|
  nodes = Array.new
  genes.map{|gene|nodes << tree.get_node_by_name(gene)}
  dist = nodes.any?(nil) ? 'nil' : tree.distance(nodes[0], nodes[1])
  distances << dist
end

puts distances.map{|i|i.to_s}.join("\t")


