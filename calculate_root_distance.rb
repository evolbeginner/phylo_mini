#! /usr/bin/env ruby


################################################
require 'bio'
require 'getoptlong'
require 'bio-nwk'


################################################
infile = nil


################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
  end
end


################################################
tree = getTreeObjs(infile).shift

tip = tree.allTips[0]

root = tree.root

puts tree.distance(root, tip)



