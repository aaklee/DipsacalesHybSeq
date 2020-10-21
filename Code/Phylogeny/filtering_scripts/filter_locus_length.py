#!/usr/bin/python

# should only be done for exon data
# supercontig data will have extraneous sequence not accounted for by the target loci (ie intron)

# input is:
#    guidance_filtering directory and
#    decimal percentage of locus length filtering threshold
# output is two lists: path of trees to cat & path of alignments to concatenate for ML analysis

import sys, os
from Bio import SeqIO
import numpy


aln_dir = os.path.abspath(sys.argv[1])
threshold = float(sys.argv[2])


target_loci = os.path.abspath('/projects/clement-lab/5-HybSeq/2-hybpiper_Dipsacales/Dipsacales_Baits_Project/Angiosperms353_targetSequences.fasta')

target_loci_lengths = {}
with open(target_loci, 'r') as inf:
    for record in SeqIO.parse(inf, 'fasta'):
        locus = record.id.split('-')[1]
        if locus not in target_loci_lengths.keys():
            target_loci_lengths[locus] = []
        target_loci_lengths[locus].append(len(record.seq))

target_loci_avg = {k: numpy.mean(v) for k, v in target_loci_lengths.items()}


drop_loci = []
keep_loci = []
for f in os.listdir(aln_dir):
    if f.endswith('.fasta'):
        locus = f[0:4]
        s = next(SeqIO.parse(os.path.join(aln_dir, f), 'fasta'))
        print(locus, len(s.seq), target_loci_avg[locus])
        if len(s.seq) < threshold*(target_loci_avg[locus]):
            drop_loci.append(locus)
        else:
            keep_loci.append(locus)

print('\n{} of {} loci < {}*(avg target locus length), {} kept\n'.format(len(drop_loci), len(target_loci_avg),  threshold, len(keep_loci)))

# generate alignment list and tree list
keep_aln = []
keep_tre = []
for f in os.listdir(aln_dir):
    if f.endswith('.fasta'):
        locus = f[0:4]
        if locus in keep_loci:
            keep_aln.append(os.path.join(aln_dir, f))
    elif 'bestTree' in f:
        locus = f.split('.')[1]
        if locus in keep_loci:
            keep_tre.append(os.path.join(aln_dir, f))

with open('{}loclen_aln_{}.txt'.format(int(threshold*100), aln_dir.split('/')[-1]), 'w+') as outf:
    for i in keep_aln:
        outf.write('{}\n'.format(i))


with open('{}loclen_tre_{}.txt'.format(int(threshold*100), aln_dir.split('/')[-1]), 'w+') as outf:
    for i in keep_tre:
        outf.write('{}\n'.format(i))


# list of filtered out loci - to avoid when filling in missing data
with open('{}loclen_filtered-out_{}.txt'.format(int(threshold*100), aln_dir.split('/')[-1]), 'w+') as outf:
    for i in drop_loci:
        outf.write('{}\n'.format(i))
