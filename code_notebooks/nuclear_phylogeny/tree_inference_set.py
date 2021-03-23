#!/usr/bin/python

# assume that {exon loci} == {supercontig loci}

# inputs:
# directories of trimmed alignments (all alignments in same top-level dir, "fasta" suffix)
# directories of trees (RAxML trees in same top-level dir, "bipartitions")
# datatype - "dna", "supercontig", etc.

# outputs:
# two text files corresponding to complete path to locus alignments/trees
# concatenated alignment, partitions, alignment summary
# concatenated trees

# next steps:
# collapse_nodes.py [optional]
# tree inference

import sys
import os

alndir = os.path.abspath(sys.argv[1])
treedir = os.path.abspath(sys.argv[2])
datatype = sys.argv[3]

paralogf = os.path.abspath('/projects/clement-lab/resources/software/DipsacalesHybSeq/code_notebooks/nuclear_phylogeny/paralogs_excluded.txt')

paralogs = []
with open(paralogf, 'r') as inf:
    paralogs = [line.strip() for line in inf]


alnf = []
treef = []
for i in os.listdir(alndir):
    if i.endswith('fasta'):
        locus = i[0:4]
        
        # skip if paralog
        if locus in paralogs:
            continue
        	
        alnf.append(os.path.join(alndir, i))

for i in os.listdir(treedir):
    if 'bipartitions.' in i:
        locus = i.split('.')[1]

        if locus in paralogs:
            continue

        treef.append(os.path.join(treedir, i))


# concatenate alignments
os.system('AMAS.py concat -f fasta -d dna -p {}.partitions.txt -t {}.gene.aln -i {}'.format(datatype, datatype, ' '.join(alnf)))
os.system('AMAS.py summary -o {}.summary.txt -f fasta -d dna -i {}.gene.aln {}'.format(datatype, datatype, ' '.join(alnf)))

# concatenate gene trees
os.system('cat {} > {}.gene.tre'.format(' '.join(treef), datatype))


# generate aln and tree file lists
with open('{}_aln_list.txt'.format(datatype), 'w+') as outf:
    outf.write('\n'.join(alnf))

with open('{}_tre_list.txt'.format(datatype), 'w+') as outf:
    outf.write('\n'.join(treef))
