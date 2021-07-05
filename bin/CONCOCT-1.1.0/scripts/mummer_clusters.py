#mummer_clusters.py - Nick Loman
# python mummer_clusters.py Contigs_15k.fasta clustering_gt1000.csv all_input_table.tsv cov_mean_sample_1122-H-Cdiff.sorted.bam cov_mean_sample_5066-H-STEC.sorted.bam > newannotations.txt

import os
import sys
import csv
from Bio import SeqIO
from collections import defaultdict

tag = sys.argv[1]
clusteringfn = sys.argv[2]
coveragefn = sys.argv[3]

##read clustering 
contigmap = {}
for ln in open(clusteringfn, "r"):
	contig, cluster_number = ln.rstrip().split(",")
	contigmap[contig] = cluster_number	

##coverage information
covmap = {}
with open(coveragefn, 'r') as csvfile:
	reader = csv.DictReader(csvfile, dialect='excel-tab')
	col1 = reader.fieldnames.index(sys.argv[4])
	col2 = reader.fieldnames.index(sys.argv[5])

with open(coveragefn, 'r') as csvfile:
	reader = csv.reader(csvfile, dialect='excel-tab')
	next(reader)
	for row in reader:
		vals = [row[0], row[1], row[2], str(sum([float(x) for x in row[col1:col2+1]]))]
		covmap[row[0]] = vals

refs = { 'NC_018658' : 'OutbreakGenome',
         'NC_018659' : 'plasmid_pESBL',
         'NC_018666' : 'plasmid_pAA',
         'NC_000913' : 'Ecoli_CoreGenome', }

contig_set = defaultdict(dict)

for chrom, label in refs.items():
	if not os.path.exists("%s_%s.delta" % (tag, chrom)):
		os.system("nucmer --prefix=%s_%s %s.fna %s" % (tag, chrom, chrom, tag))
	#os.system("show-tiling %s_%s.delta > %s_%s.tiling" % (tag, chrom, tag, chrom))
	os.system("delta-filter -i 98 %s_%s.delta > %s_%s.mummer" % (tag, chrom, tag, chrom))
	fh = open("%s_%s.mummer" % (tag, chrom,))
	while True:
		ln = fh.readline()
		if not ln: break

		if ln.startswith('>'):
			contig = ln.rstrip().split()[1]
			#contig_set[label].add(contig)
			coords = fh.readline()
			cols = coords.split(" ")
			start = cols[0]
			contig_set[label][contig] = start

print("Contig", end=' ')
for label in sorted(contig_set.keys()):
	print("\t" + label, end=' ')
print()

for rec in SeqIO.parse(tag, "fasta"):
	print("%s\t%s\t%s" % (rec.id, contigmap[rec.id], "\t".join(covmap[rec.id])), end=' ')
	for label in sorted(contig_set.keys()):
		if rec.id in list(contig_set[label].keys()):
			print("\tY\t%s" % (contig_set[label][rec.id]), end=' ')
		else:
			print("\tN\tNA", end=' ')
	print()
