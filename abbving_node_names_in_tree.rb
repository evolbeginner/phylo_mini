#! /bin/env ruby
# use abbreviation to replace the name of taxa names in a tree

##############################################################
otu_name_hash={}
p=1
tree_file = ARGV[0]


##############################################################
fh = File.open(tree_file, 'r')
while(line=fh.gets) do
  line=line.chomp
  line = line.gsub(/(?<=[\(,]) [^(:,]+/x) do |x|
    otu_name_hash[p]=1
    if x =~ /(\S+_){2,}\S*/
      $&
    else
      x =~ /^(\S)\S+_(.{,3}).+/
      p=$1+$2
      if otu_name_hash.include?($1+$2)
        $&
      else
        $1+$2
      end
    end
  end
  puts line
end


