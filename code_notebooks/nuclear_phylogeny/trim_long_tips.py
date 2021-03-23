#!/usr/bin/python

import sys
import os
import dendropy
from Bio import SeqIO
import argparse
#from matplotlib import pyplot
import numpy

OUTGROUP_TAXA = ['Ilex_crenata', 'Paracryphia_alticola'] # edit this as necessary

def trim(tre, cutoff):
    cutoff_sequences = []
    edge_lens = []
    #tre.print_plot()
    for node in tre:
        if node.num_child_nodes() == 0:
            if node.edge_length > cutoff:
                cutoff_sequences.append(str(node.taxon))
            edge_lens.append(float(node.edge_length))

        # trim tip branches if they're much longer than their neighbors?
        #elif node.num_child_nodes() == 2:
        #    print(node.child_nodes())

    cutoff_sequences = [i.replace(' ', '_') for i in cutoff_sequences]
    cutoff_sequences = [i.replace('\'', '') for i in cutoff_sequences]

    return edge_lens, cutoff_sequences


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tree_dir', help='path to directory with trees')
    parser.add_argument('-a', '--aln_dir', help='path to directory with alignments for trees')
    parser.add_argument('-c', '--cutoff', type=float, help='tip branch length cutoff, e.g. 1.0' )
    parser.add_argument('-p', '--len_plot', action='store_true', help='include if you want to generate a plot of tip branch lengths')
    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    tree_dir = os.path.abspath(args.tree_dir)
    aln_dir = os.path.abspath(args.aln_dir)
    cutoff = float(args.cutoff)
    len_plot = args.len_plot


    # trim long terminal branches
    print('\n~~ TRIMMING ~~')

    # obtain sequences that should be trimmed from dataset
    to_trim = {}
    tip_edge_lens = {}
    for i in os.listdir(tree_dir):
        if 'bestTree' in i:
            locus = i.split('.')[1]
            tre = dendropy.Tree.get(path=os.path.join(tree_dir, i), schema='newick')
            tip_edge_lens[locus], to_trim[locus] = trim(tre, cutoff)


    for i in os.listdir(aln_dir):
        if i.endswith('.fasta'):
            locus = i[0:4]

            # filter out trimmed sequences from original fasta file
            keep_seqs = []
            seqct = 0
            cutct = 0


            with open(os.path.join(aln_dir, i), 'r') as inf:
                for record in SeqIO.parse(inf, 'fasta'):

                    # don't trim outgroup taxa
                    if record.id in OUTGROUP_TAXA:
                        '{} (outgroup) has a long branch'.format(record.id)
                        pass
                    elif locus in to_trim.keys() and record.id in to_trim[locus]:
                        cutct += 1
                        continue
                    keep_seqs.append(record)
                    seqct += 1

            if locus in to_trim.keys():
                print('{}\t{}->{}\t-{}\t{}'.format(locus, seqct, seqct-cutct, len(to_trim[locus]), ','.join(to_trim[locus])))
            else:
                print('{}\t{}->{}\tNONE'.format(locus, seqct, seqct-cutct))

            # write (sequences minus trimmed sequences) to output fasta file
            newfn = '.'.join([i.split('.')[0], 'trimmed', 'fasta'])
            with open(newfn, 'w+') as outf:
                SeqIO.write(keep_seqs, outf, 'fasta')

    print('~~ TRIMMING COMPLETE ~~\n')



    # calculate stats
    loci = list(tip_edge_lens.keys())
    lens_list = list(tip_edge_lens.values())

    lens = []
    for i in lens_list:
        for j in i:
            lens.append(j)

    mean = numpy.mean(lens)
    median = numpy.median(lens)
    stdev = numpy.std(lens)
    min = numpy.amin(lens)
    max = numpy.amax(lens)

    print('\n~~ TERMINAL BRANCH LEN STATS ~~')
    print('mean\t=\t{} +/- {}'.format(mean, stdev))
    print('upper\t=\t{}'.format(mean+stdev))
    print('lower\t=\t{}'.format(mean-stdev))
    print('median\t=\t{}'.format(median))
    print('min\t=\t{}'.format(min))
    print('max\t=\t{}'.format(max))


    # this doesn't work
    #if len_plot:
    #    print('generating plot')
    #    fig, ax = pyplot.subplots()
    #    for i in range(len(loci)):
    #        ax.scatter(loci, lens)
    #    fig.savefig('plot.pdf')



if __name__ == '__main__':
    main()
