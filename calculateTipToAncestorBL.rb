#! /usr/bin/env ruby


####################################################
dir = File.dirname($0)
$: << File.join(dir, "treeRubyLib")


####################################################
require 'getoptlong'
require 'bio'

require 'tree'


####################################################
infile = nil


####################################################
if __FILE__ == $0
  opts = GetoptLong.new(
    ['-i', GetoptLong::REQUIRED_ARGUMENT],
  )


  opts.each do |opt, value|
    case opt
      when /^-i$/
        infile = value
    end
  end


  ####################################################
  trees = getTreeObjs(infile)
  trees.each do |tree|
    tips = tree.nodes.select{|node|node.isTip?(tree)}
    ancestor_node = tips.inject{|ancestor, tip| tree.lowest_common_ancestor(tip, ancestor)}
    begin
      puts tips.map{|tip| [tip.name, tree.distance(tip, ancestor_node)].join("\t")}.join("\n")
    rescue Exception => e
      puts "#{e}"
    end
  end

end

