#! /bin/env ruby


##############################################################
require 'getoptlong'

require 'util'


#############################################################
infile = nil
must_file = nil
prior_file = nil
is_only = false

must = Hash.new
prior = Hash.new
clstrs = Hash.new{|h1,k1|h1[k1]=Hash.new{|h2,k2|h2[k2]=[]}}


#############################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--must', GetoptLong::REQUIRED_ARGUMENT],
  ['--prior', GetoptLong::REQUIRED_ARGUMENT],
  ['--only', GetoptLong::REQUIRED_ARGUMENT],
  ['--exclude_prior', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '--must'
      must_file = value
    when '--prior'
      prior_file = value
    when '--exclude_prior'
      exclude_prior_file = value
    when '--only'
      is_only = true
      prior_file = value
  end
end


#############################################################
must = read_list(must_file).keys unless must_file.nil?
prior = read_list(prior_file).keys unless prior_file.nil?
#exclude_prior = read_list(exclude_prior_file).keys unless exclude_prior_file.nil?

in_fh = infile == '-' ? STDIN : File.open(infile, 'r')


#############################################################
in_fh.each_line do |line|
  line.chomp!
  species, clstr = line.split("\t")
  clstr = clstr.to_i
  if clstr == -1
    s = clstrs[:independent].size
    clstrs[:independent][s-1] << species
  else
    clstrs[:others][clstr] << species
  end
end

in_fh.close if infile != '-'


#############################################################
if is_only
  clstrs[:others].each_pair do |clstr, species_names|
    if species_names.any?{|i| prior.include?(i) }
      puts species_names.select{|i| prior.include?(i) }[0]
    end
  end

  clstrs[:independent].each_pair do |clstr, species_names|
    if species_names.any?{|i| prior.include?(i) }
      puts species_names.select{|i| prior.include?(i) }[0]
    end
  end

  exit
end


clstrs.each_pair do |type, h|
  h.each_pair do |clstr, species_names|
    if type == :others
      if species_names.any?{|i| must.include?(i) }
        puts species_names.select{|i| must.include?(i) }.join("\n")
      elsif species_names.any?{|i| prior.include?(i) }
        puts species_names.select{|i| prior.include?(i) }[0]
      elsif species_names.any?{|i| i !~ /^proteo/}
        puts species_names.select{|i|i !~ /^proteo/}[0]
      else
        puts species_names[0]
        ;
      end
    elsif type == :independent
      puts species_names[0]
    end
  end
end


