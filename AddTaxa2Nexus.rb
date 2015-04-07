#! /bin/env ruby

require 'getoptlong'

###################################################
nexus_file=nil
taxa=nil

opts=GetoptLong.new(
  ['-i', '--in', GetoptLong::REQUIRED_ARGUMENT],
  ['-t', '--taxa', GetoptLong::REQUIRED_ARGUMENT],
)

opts.each do |opt,value|
  case opt
    when '-i', '--in'
      nexus_file=value
    when '-t', '--taxa'
      taxa=value
  end
end

###################################################
fh=File.open(nexus_file,'r')
while(line=fh.gets) do
  line.chomp!
  puts line
  if line =~ /^BEGIN TREES;$/ then
    while(line=fh.gets) do
      line.chomp!
      if line =~ /^\s+\d+\s+/ then
        line.sub!(/^(\s+\d+\s+)(.+)/) do
          $1+[taxa,$2].join("_")
        end
      end
      puts line
    end
  end
end
fh.close

=begin
BEGIN TREES;
      TITLE 'Trees from "original.above50.nexus"';
      LINK TAXA = Taxa;
         TRANSLATE
            1    Arabidopsis_thaliana_AT3G59590,
            2    Arabidopsis_thaliana_AT3G59610,
=end

