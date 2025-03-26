#! /usr/bin/env ruby


##################################################
require 'getoptlong'
require 'Dir'


##################################################
indir = nil
outdir = nil
is_force = false


##################################################
opts = GetoptLong.new(
  ['--indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT]
)

opts.each do |opt, value|
  case opt
    when '--indir'
      indir = value
    when '--outdir'
      outdir = value
    when '--force'
      is_force = true
  end
end


##################################################
mkdir_with_force(outdir, is_force)

tree_outdir = File.join(outdir, "phylo", "iqtree")
mkdir_with_force(tree_outdir, is_force)

# combined.aln
aln_outfile = File.join(outdir, "combined.aln")
`cp #{indir}/combined.aln #{aln_outfile}`

# treefile
tree_outfile = File.join(tree_outdir, "REV.treefile")
`cp #{indir}/phylo/iqtree/iqtree.treefile #{tree_outfile}`
`nw_topology -Ib #{tree_outfile} | sponge #{tree_outfile}`

# best_scheme
scheme_infile = File.join(indir, "phylo/iqtree/iqtree.best_scheme.nex")
if File.exists?(scheme_infile)
  scheme_outfile = File.join(tree_outdir, "REV.best_scheme.nex")
  `cp #{scheme_infile} #{scheme_outfile}`
else
  puts "scheme file does not exist"
end


