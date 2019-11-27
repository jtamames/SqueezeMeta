#!/usr/bin/env ruby

# The MIT License (MIT)
# Copyright (c) 2016 Alexander J Probst

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# https://github.com/AJProbst/sngl_cp_gn

#d = ARGV[0]

input_file = ARGV[1]
output_dir = File.dirname(input_file)

datab_in = ARGV[2]
datab_out = "#{File.dirname ARGV[1]}" + "/" + "#{File.basename ARGV[2]}"
puts "database name of all proteins is #{datab_in}"

db_name = ARGV[3]
puts "database name of SCGs is #{db_name}"

db_lookup = ARGV[4]
puts "database lookup is #{db_lookup}"

threads = ARGV[5]

#build databases
full_db = system "makeblastdb -in #{datab_in} -dbtype prot -out #{datab_out}"
abort "makeblastdb did not work for #{datab_in}, please check your input file" unless full_db

# find SCG candidates
puts "finding SCG candidates..."
input_blast_database = system "makeblastdb -in #{input_file} -dbtype prot"
input_blast_out = File.join(output_dir,File.basename(input_file) + ".findSCG.b6")
abort "makeblastdb did not work for #{input_file}, please check your input file" unless input_blast_database
input_blast_ok = system "blastp -db #{input_file} -query #{db_name} -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' -evalue 0.01 -out #{input_blast_out} -num_threads #{threads}"
system "rm #{input_file}.psq #{input_file}.pin #{input_file}.phr"
abort "blast did not work, please check your input file." unless input_blast_ok

input_blast_out_whitelist = File.join(output_dir,File.basename(input_file) + ".findSCG.b6.whitelist")
system "awk '{print$2}' #{input_blast_out} | sort -u > #{input_blast_out_whitelist}"
scg_candidates = File.join(output_dir,File.basename(input_file) + ".scg.candidates.faa")
system "pullseq -i #{input_file} -n #{input_blast_out_whitelist} > #{scg_candidates}"
system "rm #{input_blast_out_whitelist}"

# verify SCGs by blasting against all proteins of all genomes
puts "verifying selected SCGs..."
db_blast_out = File.join(output_dir,File.basename(input_file) + ".all.b6")
db_blast_ok = system "blastp -db #{datab_out} -query #{scg_candidates} -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' -evalue 0.00001 -out #{db_blast_out} -max_target_seqs 1 -num_threads #{threads}"
abort "verifying blast did not work" unless db_blast_ok
system "rm #{datab_out}.psq #{datab_out}.pin #{datab_out}.phr"
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


