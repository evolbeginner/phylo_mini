#! /bin/env ruby


###############################################################
require 'getoptlong'
require 'bio-nwk'


###############################################################
def get_ordered_tips(infile)
  ordered_tips = Array.new
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line.scan(/([^(),: ]+):/) do |i|
      ordered_tips << i[0]
    end
  end
  in_fh.close
  return(ordered_tips)
end


def get_child(tree, node, bls)
  if node != tree.root
    puts tree.distance(node, tree.parent(node))
    bls << tree.distance(node, tree.parent(node))
  end

  if node.isTip?(tree)
    #p node.name
    ;
  else
    children = tree.children(node)
    ordered_children = children.sort_by{|child| tree.tips(child).map{|i|ORDERED_TIPS.index(i.name)}.min}
    #p [ordered_children, node]
    ordered_children.each do |i|
      get_child(tree, i, bls)
    end
  end
end


###############################################################
infile = nil
ref_tree_file = nil

ordered_tips = Array.new


###############################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--ref', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '--ref'
      ref_tree_file = value
  end
end


###############################################################
ORDERED_TIPS = get_ordered_tips(ref_tree_file)

trees = getTreeObjs(infile)


###############################################################
trees.each do |tree|
  bls = Array.new
  root = tree.root
  get_child(tree, root, bls)
  p bls
end


