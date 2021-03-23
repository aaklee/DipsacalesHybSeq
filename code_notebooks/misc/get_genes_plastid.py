#!/usr/bin/python

import sys
import os
import pandas
from Bio import SeqIO
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--list', help='list of alignments used in final tree', required=True)
    parser.add_argument('-a', '--abs_taxa', help='OPTION 1: how many taxa you want to include (do not use with -p, -n)')
    parser.add_argument('-p', '--percent_loci', help='OPTION 2: percent of total loci to include (do not use with -n, -a)')
    parser.add_argument('-n', '--num_loci', help='OPTION 3: how many loci you want to include (do not use with -p, -a)')
    parser.add_argument('-s', '--sum_bp', help='OPTION 4: how many bp you want to include (do not use with others)')

    args = parser.parse_args()

    if args.abs_taxa:
        abs_taxa = int(args.abs_taxa)
    else:
        abs_taxa = False

    if args.percent_loci:
        percent_loci = float(args.percent_loci)
    else:
        percent_loci = False

    if args.num_loci:
        num_loci = int(args.num_loci)
    else:
        num_loci = False

    if args.sum_bp:
        sum_bp = int(args.sum_bp)
    else:
        sum_bp = False


    aln = []
    with open(os.path.abspath(args.list), 'r') as inf:
        for line in inf:
            aln.append(line.strip())

    aln_dict = {}
    for a in aln:
        locus = a.split('/')[-1][0:4]

        seq_list = []
        with open(a, 'r') as inf:
            seq_list = list(SeqIO.parse(inf, 'fasta'))

        if len(seq_list) <= 3:
            continue

        for i in range(len(seq_list)):
            seq_list[i].description = ''

        num_seqs = len(seq_list)
        aln_len = len(seq_list[0].seq)

        aln_dict[locus] = [a, num_seqs, aln_len]

    aln_df = pandas.DataFrame.from_dict(aln_dict, orient='index', columns=['path', 'num_seqs', 'aln_len'])
    aln_df = aln_df.sort_values(by=['num_seqs'], ascending=False)

    if abs_taxa:
        aln_df_filt = aln_df[aln_df['num_seqs'] >= abs_taxa]
        total_len = aln_df_filt['aln_len'].sum()
        print(aln_df_filt[['num_seqs', 'aln_len']])

        print(aln_df_filt.shape)

        print('TOTAL LOCI:', len(aln_df_filt))
        print('TOTAL LEN:', total_len)
        print('MEAN LEN:', total_len / len(aln_df_filt))

        alns = aln_df_filt['path']
        aln_cat = [i for i in alns]
        print(len(aln_cat))

        os.system('mkdir {}taxa'.format(abs_taxa))
        os.system('AMAS.py concat -i {} -f fasta -d dna -t concatenated.fasta'.format(' '.join(aln_cat)))
        os.system('cp {} ./'.format(' '.join(aln_cat)))
        os.system('mv *.fasta *.txt {}taxa'.format(abs_taxa))


    elif percent_loci:
        num = 60 * percent_loci
        perc = int(percent_loci * 100)

        aln_df_filt = aln_df[aln_df['num_seqs'] >= num]
        total_len = aln_df_filt['aln_len'].sum()
        print(aln_df_filt[['num_seqs', 'aln_len']])
        print('TOTAL LOCI:', len(aln_df_filt))
        print('TOTAL LEN:', total_len)

        alns = aln_df_filt['path']
        aln_cat = [i for i in alns]

        os.system('mkdir top_{}perc'.format(perc))
        os.system('AMAS.py concat -i {} -f fasta -d dna -t concatenated.fasta'.format(' '.join(aln_cat)))
        os.system('cp {} ./'.format(' '.join(aln_cat)))
        os.system('mv *.fasta *.txt top_{}perc'.format(perc))

    elif num_loci:
        aln_df_filt = aln_df.head(num_loci)
        total_len = aln_df_filt['aln_len'].sum()
        print(aln_df_filt[['num_seqs', 'aln_len']])
        print('TOTAL LEN:', total_len)
        print('TOTAL LOCI:', len(aln_df_filt))

        alns = aln_df_filt['path']
        aln_cat = [i for i in alns]

        os.system('mkdir top_{}'.format(num_loci))
        os.system('AMAS.py concat -i {} -f fasta -d dna -t concatenated.fasta'.format(' '.join(aln_cat)))
        os.system('cp {} ./'.format(' '.join(aln_cat)))
        os.system('mv *.fasta *.txt top_{}'.format(num_loci))

    elif sum_bp:
        sum = 0
        alns = []
        for i, j in aln_df.iterrows():
            sum += j['aln_len']
            alns.append(j['path'])
            if sum >= sum_bp:
                break
        total_len = sum
        total_loci = len(alns)
        print(aln_df[['num_seqs', 'aln_len']])
        print('FINAL LEN:', total_len)
        print('FINAL LOCI:', total_loci)

        aln_cat = [i for i in alns]

        print(' '.join(aln_cat))

        os.system('mkdir {}bp'.format(sum))
        os.system('AMAS.py concat -i {} -f fasta -d dna -t concatenated.fasta'.format(' '.join(aln_cat)))
        os.system('cp {} ./'.format(' '.join(aln_cat)))
        os.system('mv *.fasta *.txt {}bp'.format(sum))



if __name__ == '__main__':
    main()
