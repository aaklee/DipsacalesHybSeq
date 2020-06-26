#!/usr/bin/python

# count total loci, total taxa from a directory containing fasta alignments

import sys, os
from Bio import SeqIO

alndir = os.path.abspath(sys.argv[1])

loci = []
taxa = []
for f in os.listdir(alndir):
    if f.endswith('.FNA') or f.endswith('.fasta'):
        locus = f[0:4]
        loci.append(locus)
        with open(os.path.join(alndir, f), 'r') as inf:
            for record in SeqIO.parse(inf, 'fasta'):
                taxa.append(record.id[0:10])

print('num loci: {}'.format(len(set(loci))))
print('num taxa: {}'.format(len(set(taxa))))
