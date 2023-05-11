#! /bin/env ruby


############################################################
require 'getoptlong'
require 'parallel'
require 'colorize'

require 'util'
require 'Dir'


############################################################
COLOR_TIP = File.expand_path("~/tools/self_bao_cun/phylo_mini/graph/color_tip.R")

indir = nil
cpu = 1
font_size = 0.8
color_file = nil
layout = 'circular'
is_bl = true
type = nil
suffix = nil
add_argu = ''


############################################################
opts = GetoptLong.new(
  ['--indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
  ['--layout', GetoptLong::REQUIRED_ARGUMENT],
  ['-s', GetoptLong::REQUIRED_ARGUMENT],
  ['-a', GetoptLong::REQUIRED_ARGUMENT],
  ['--type', GetoptLong::REQUIRED_ARGUMENT],
  ['--add_argu', GetoptLong::REQUIRED_ARGUMENT],
  ['--no_bl', GetoptLong::NO_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '--indir'
      indir = value
    when '--cpu'
      cpu = value.to_i
    when '--layout'
      layout = value
    when '-s'
      font_size = value.to_f
    when '-a'
      color_file = value
    when '--no_bl'
      is_bl = false
    when '--add_argu'
      add_argu = value
    when '--type'
      type = value
  end
end


############################################################
case type
  when /iqtree/
    suffix = 'treefile'
  when /fasttree/i
    suffix = 'fasttree.tre'
end


if type.nil?
  STDERR.puts "--type has to be specified!".colorize(:red)
  exit 1
end


############################################################
infiles = read_infiles(indir)
Parallel.map(infiles, in_processes: cpu) do |infile|
  c = getCorename(infile)
  puts c
  treefile = File.join(infile, c + '.' + suffix)
  outfile = File.join(infile, [c, type, layout, 'pdf'].join('.'))
  bl_argu = "--no_bl" if not is_bl
  p "Rscript #{COLOR_TIP} -a #{color_file} -t #{treefile} -o #{outfile} #{bl_argu} --layout #{layout} -s #{font_size} #{add_argu}"
  `Rscript #{COLOR_TIP} -a #{color_file} -t #{treefile} -o #{outfile} #{bl_argu} --layout #{layout} -s #{font_size} #{add_argu}`
end


