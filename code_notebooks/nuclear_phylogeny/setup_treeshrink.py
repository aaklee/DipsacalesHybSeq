#!/usr/bin/python

# set up treeshrink for a per-species analysis using alignments
# this helps to just filter out any tips that were trimmed by treeshrink

# here is the next step:
# python /projects/clement-lab/resources/software/TreeShrink-1.3.8b/run_treeshrink.py -i dna -t input.tre -a input.fasta -f -b 0.1 > dna.log

import sys
import os

alns = os.path.abspath(sys.argv[1])
trees = os.path.abspath(sys.argv[2])
outdir = sys.argv[3]

loci = {}
for i in os.listdir(alns):
	if i.endswith('fasta'):
		locus = i.split('_')[0]
		loci[locus] = [os.path.join(alns, i)]

for i in os.listdir(trees):
	if 'bipartitions.' in i:
		locus = i.split('.')[1]
		loci[locus].append(os.path.join(trees, i))

os.system('mkdir {}'.format(outdir))
for k,v in loci.items():
	os.system('mkdir {}/{}'.format(outdir, k))
	os.system('cp {} {}/{}/input.fasta'.format(v[0], outdir, k))
	os.system('cp {} {}/{}/input.tre'.format(v[1], outdir, k))


