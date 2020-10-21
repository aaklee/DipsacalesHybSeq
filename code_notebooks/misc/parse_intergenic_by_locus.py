#!/usr/bin/python

import sys
import os
from Bio import SeqIO
import pandas

def main():
    min_taxa = int(sys.argv[1])

    fastas = "/projects/clement-lab/5-HybSeq/hybpiper_plastids/partition_extractions/intergenic/extracted_regions"
    tallies = "/projects/clement-lab/5-HybSeq/hybpiper_plastids/partition_extractions/intergenic/tallied_regions.txt"

    df = pandas.read_csv(tallies, sep='\t', header=0)
    df_sub = df.loc[df['num_species'] >= min_taxa]

    include_loci = {i: [] for i in df_sub['intergenic_region'].values.tolist()}

    for i in os.listdir(fastas):
        if i.endswith('fasta'):
            for record in SeqIO.parse(os.path.join(fastas, i), 'fasta'):
                region = record.id.split('.')[1]
                if region in include_loci.keys():
                    include_loci[region].append(record)

    for i, j in include_loci.items():
        newf = '{}.fasta'.format(i)
        SeqIO.write(j, newf, 'fasta')

if __name__ == '__main__':
    main()
