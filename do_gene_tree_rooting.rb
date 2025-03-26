#! /usr/bin/env ruby


##############################################
require 'getoptlong'
require 'parallel'
require 'find'


##############################################
MAD = File.expand_path("~/software/phylo/MAD/mad.py")
LSD = File.expand_path("~/software/phylo/lsd-0.2/src/lsd")

ACTIVATE='conda activate py3.7; export PYTHONPATH=""'


##############################################
indir = nil
cpu = 4
suffixes = %w[tre treefile]

treefiles = Array.new


##############################################
opts = GetoptLong.new(
  ['--indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
  ['--suffix', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '--indir'
      indir = value
    when '--cpu'
      cpu = value.to_i
    when '--suffix'
      suffixes = value.split(',')
  end
end


if indir.nil?
	STDERR.puts "indir has to be given!"
	exit 1
end


##############################################
Find.find(indir) do |path|
  treefiles << path if suffixes.any?{|suffix| path =~ /\.#{suffix}$/}
end


Parallel.map(treefiles, in_threads: cpu) do |treefile|
  puts treefile
	`python #{MAD} -t -n #{treefile}`
  `#{LSD} -r a -i #{treefile} -o #{treefile}.lsd`
  #p "#{LSD} -r a -i #{treefile} -o #{treefile}.lsd"
end


