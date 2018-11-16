#!/usr/bin/env ruby -w

# The MIT License (MIT)
# Copyright (c) 2016 Alexander J Probst

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# https://github.com/AJProbst/sngl_cp_gn

#1: $search_engine name
#2: $proteins
#3: $DIR\/db/bac.all.faa
#4: $DIR\/db/bac.scg.faa
#5: $DIR\/db/bac.scg.lookup
#6: $threads

d = ARGV[0]

input_file = ARGV[1]
output_dir = File.dirname(input_file)

datab = ARGV[2]
db_all = File.dirname(input_file) + "/all_prot"
puts "database name of all proteins is #{datab}"

db_name = ARGV[3]
puts "database name of SCGs is #{db_name}"

db_lookup = ARGV[4]
puts "database lookup is #{db_lookup}"

threads = ARGV[5]

#build databases
full_db = system "#{d} makedb --in #{datab} -d #{db_all}.dmnd"
abort "makeblastdb did not work for #{datab}, please check your input file" unless full_db

# find SCG candidates
puts "finding SCG candidates..."
input_blast_database = system "#{d} makedb --in #{input_file} -d #{input_file}.dmnd"
input_blast_out = File.join(output_dir,File.basename(input_file) + ".findSCG.b6")
abort "makeblastdb did not work for #{input_file}, please check your input file" unless input_blast_database
input_blast_ok = system "#{d} blastp --query #{db_name} --db #{input_file}.dmnd --max-target-seqs 0 --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore --out #{input_blast_out} --evalue 0.01 --threads #{threads}"
system "rm #{input_file}.dmnd"
abort "blast did not work, please check your input file." unless input_blast_ok

input_blast_out_whitelist = File.join(output_dir,File.basename(input_file) + ".findSCG.b6.whitelist")
system "awk '{print$2}' #{input_blast_out} | sort -u > #{input_blast_out_whitelist}"
scg_candidates = File.join(output_dir,File.basename(input_file) + ".scg.candidates.faa")
system "pullseq -i #{input_file} -n #{input_blast_out_whitelist} > #{scg_candidates}"
system "rm #{input_blast_out_whitelist}"

# verify SCGs by blasting against all proteins of all genomes
puts "verifying selected SCGs..."
db_blast_out = File.join(output_dir,File.basename(input_file) + ".all.b6")
db_blast_ok = system "#{d} blastp --query #{scg_candidates} --db #{db_all} --evalue 0.00001 --threads #{threads} --out #{db_blast_out} --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore --max-target-seqs 1"
abort "verifying blast did not work" unless db_blast_ok
system "rm #{db_all}.dmnd"
puts "starting annotations of single copy cogs..."

# Read db_lookup
lookup_h = {}
File.open(db_lookup).each do |line|
  sbj, annotation = line.chomp.split
  lookup_h[sbj]=annotation
end

# now compare and print
File.open(File.join(output_dir,File.basename(input_file)+".scg"), "w") do |file|
  File.open(db_blast_out).each do |line|
    next if line =~ /^#/
    line.chomp!
    temp = line.split(/\t/)
    query, sbjct = temp[0], temp[1]
    aln_len, sbjct_len = temp[3], temp[5]
    if lookup_h[sbjct] && aln_len > (sbjct_len*0.5)
      file.puts "#{query.split[0]}\t#{lookup_h[sbjct]}"
    end
  end
end

puts "successfully finished"
