#!/usr/bin/python

import sys
import os
from Bio import SeqIO
import pandas

def main():
    dir = os.path.abspath(sys.argv[1])

    intergenic = {}
    for i in os.listdir(dir):
        if i.endswith('fasta'):
            species = '_'.join(i.split('_')[:2])
            for record in SeqIO.parse(os.path.join(dir, i), 'fasta'):
                region = record.id.split('.')[1]

                if region not in intergenic:
                    intergenic[region] = []

                intergenic[region].append(species)

    with open('tallied_regions.txt', 'w+') as outf:
        outf.write('intergenic_region\tnum_species\tspecies\n')
        for i, j in intergenic.items():
            outf.write('{}\t{}\t{}\n'.format(i, len(j), ','.join(j)))



if __name__ == '__main__':
    main()
