#!/usr/bin/python

import sys
import os
from Bio import SeqIO

def main():
    aln_dir = os.path.abspath(sys.argv[1]) # path to PASTA alignments
    cutoff = float(sys.argv[2]) # phyutility -clean cutoff

    for i in os.listdir(aln_dir):
        if os.path.isdir(os.path.join(aln_dir,i)):
            locus = i

            for j in os.listdir(os.path.join(aln_dir, i)):
                if j.endswith('.aln'):
                    c = int(cutoff * 100)
                    os.system('java -jar /projects/clement-lab/resources/software/phyutility/phyutility.jar -clean {} -in {} -out {}_{}.fasta'.format(cutoff, os.path.join(aln_dir,i,j), i, c))


    # remove completely undetermined sequences, alignments with <4 taxa
    for i in os.listdir(os.getcwd()):
        if i.endswith('.fasta'):
            seqs = []
            for record in SeqIO.parse(i, 'fasta'):
                basect = 0
                for base in record.seq:
                    if base != '-':
                        break
                    basect += 1

                if basect == len(record.seq):
                    print('removing {} from {}'.format(record.id, i))
                else:
                    seqs.append(record)

            if len(seqs) < 4:
                print('< 4 seqs, removing {}'.format(i))
                os.system('rm {}'.format(i))
            else:
                with open(i, 'w+') as outf:
                    SeqIO.write(seqs, outf, 'fasta')






if __name__ == '__main__':
    main()
