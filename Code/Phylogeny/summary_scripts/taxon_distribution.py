#!/usr/bin/python

import sys
import os
from Bio import SeqIO
import pandas


aln_dir = os.path.abspath(sys.argv[1])
tre_dir = os.path.abspath(sys.argv[2])

aln_all = []
tre_loci = []
for f in os.listdir(aln_dir):
    if f.endswith('.fasta') or f.endswith('.FNA'):
        aln_all.append(os.path.join(aln_dir, f))

for f in os.listdir(tre_dir):
    if 'bestTree' in f:
        locus = f.split('.')[1]
        tre_loci.append(locus)

aln_f = [x for x in aln_all if x.split('/')[-1][0:4] in tre_loci]


# distribution of taxa among loci
taxon_distribution = {}
for f in aln_f:
    with open(f, 'r') as inf:
        locus = f.split('/')[-1][0:4]
        taxon_distribution[locus] = []
        for record in SeqIO.parse(f, 'fasta'):
            taxon_distribution[locus].append(record.id)

with open('taxon_distribution.txt', 'w+') as outf:
    outf.write('locus\ttaxa\n')
    for k, v in taxon_distribution.items():
        outf.write('{}\t{}\n'.format(k, ','.join(v)))


# distribution of loci among taxa
locus_distribution = {}
for k, v in taxon_distribution.items():
    for s in v:
        if s not in locus_distribution.keys():
            locus_distribution[s] = []
        locus_distribution[s].append(k)

with open('locus_distribution.txt', 'w+') as outf:
    outf.write('taxon\tloci\n')
    for k, v in locus_distribution.items():
        outf.write('{}\t{}\n'.format(k, ','.join(v)))


# binary matrix?
#td_df = pandas.DataFrame(taxon_distribution)
#print(td_df)
