#!/usr/bin/python

import sys
import os
import argparse
import pandas
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-sa', '--supercontig_aln', help='path to directory with supercontig alignments')
    parser.add_argument('-st', '--supercontig_tre', help='path to directory with supercontig trees')
    parser.add_argument('-ea', '--exon_aln', help='path to directory with exon alignments')
    parser.add_argument('-et', '--exon_tre', help='path to directory with exon trees')
    parser.add_argument('-a', '--avoid_loci', help='output "filtered_loci" from filtering by locus length')
    parser.add_argument('-b', '--bootstrap', action='store_true', help='include if you want to generate concatenated bootstrap trees/accompanying lists' )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    supercontig_aln = os.path.abspath(args.supercontig_aln)
    supercontig_tre = os.path.abspath(args.supercontig_tre)
    exon_aln = os.path.abspath(args.exon_aln)
    exon_tre = os.path.abspath(args.exon_tre)


    # these are loci for which the recovered exons did reach a certain percentage of the target locus
    # if found in supercontig dataset, flagged but included due to flanking regions providing additional information
    # if found in exon dataset, ignored.
    avoid = []
    if args.avoid_loci:
        with open(args.avoid_loci, 'r') as inf:
            for line in inf:
                avoid.append(line.strip())

    flagged_loci = []
    ignored_loci = []



    # get which loci are inverse between the supercontig/exon trees and alignments
    # (some loci will have dropped out due to RAxML only taking loci with >=4 taxa)
    # ensures consistency between coalescent and concatenated analyses
    final_loci = []
    for f in os.listdir(supercontig_tre):
        if args.bootstrap:
            if 'bipartitions.' in f:
                locus = f.split('.')[1]
                if locus not in final_loci:
                    final_loci.append(locus)

        else:
            if 'bestTree' in f:
                locus = f.split('.')[1]
                if locus not in final_loci:
                    final_loci.append(locus)

    for f in os.listdir(exon_tre):
        if args.bootstrap:
            if 'bipartitions.' in f:
                locus = f.split('.')[1]
                if locus not in final_loci:
                    final_loci.append(locus)

        else:
            if 'bestTree' in f:
                locus = f.split('.')[1]
                if locus not in final_loci:
                    final_loci.append(locus)

    print(len(final_loci))


    # first take all supercontigs (supercontig-only)
    # second take all exons for which there were no supercontig dataset (supercontig+exon)
    # third take all exons (exon-only)

    to_cat = {'supercontig':[], 'inverse':[], 'exon':[]} # tre files to concatenate
    to_aln = {'supercontig':[], 'inverse':[], 'exon':[]} # fasta files to concatenate
    loci = {'supercontig':[], 'inverse':[], 'exon':[]} # which loci are being included?
    # NOTE: taxa are not being accounted for yet -- should be the same between coalescent/concat analyses, so we'll get those later

    for f in os.listdir(supercontig_tre):
        if args.bootstrap:
            if 'bipartitions.' in f:
                locus = f.split('.')[1]

                if locus not in avoid:
                    if locus in final_loci: # ignores a locus if it should be filtered out based on locus length
                        to_cat['supercontig'].append(os.path.join(supercontig_tre, f))
                        loci['supercontig'].append(locus)
                else:
                    flagged_loci.append(locus)

        else:
            if 'bestTree' in f:
                locus = f.split('.')[1]

                if locus not in avoid:
                    if locus in final_loci: # ignores a locus if it should be filtered out based on locus length
                        to_cat['supercontig'].append(os.path.join(supercontig_tre, f))
                        loci['supercontig'].append(locus)
                else:
                    flagged_loci.append(locus)

    for f in os.listdir(supercontig_aln):
        if f.endswith('.fasta'):
            locus = f[0:4]

            if locus in loci['supercontig']:
                if locus not in avoid:
                    if locus in final_loci:
                        to_aln['supercontig'].append(os.path.join(supercontig_aln, f))


    for f in os.listdir(exon_tre):
        if args.bootstrap:
            if 'bipartitions.' in f:
                locus = f.split('.')[1]

                if locus not in avoid:
                    if locus in final_loci:
                        if locus not in loci['supercontig']: # exon locus not yet accounted for by supercontig data
                            to_cat['inverse'].append(os.path.join(exon_tre, f))
                            loci['inverse'].append(locus)
                        to_cat['exon'].append(os.path.join(exon_tre, f))
                        loci['exon'].append(locus)
                else:
                    ignored_loci.append(locus)

        else:
            if 'bestTree' in f:
                locus = f.split('.')[1]

                if locus not in avoid:
                    if locus in final_loci:
                        if locus not in loci['supercontig']: # exon locus not yet accounted for by supercontig data
                            to_cat['inverse'].append(os.path.join(exon_tre, f))
                            loci['inverse'].append(locus)
                        to_cat['exon'].append(os.path.join(exon_tre, f))
                        loci['exon'].append(locus)
                else:
                    ignored_loci.append(locus)

    for f in os.listdir(exon_aln):
        if f.endswith('.fasta'):
            locus = f[0:4]
            if locus not in avoid and locus in final_loci:
                if locus not in loci['supercontig']:
                    to_aln['inverse'].append(os.path.join(exon_aln, f))
                to_aln['exon'].append(os.path.join(exon_aln, f))

    # generate files of trees for coalescent analysis
    print('\n~ CONCATENATE GENE TREES ~')
    os.system('cat {} > supercontig.gene.tre'.format(' '.join(to_cat['supercontig'])))
    os.system('cat {} {} > supercontig-exon.gene.tre'.format(' '.join(to_cat['supercontig']), ' '.join(to_cat['inverse'])))
    os.system('cat {} > exon.gene.tre'.format(' '.join(to_cat['exon'])))
    print('~ CONCATENATE GENE TREES COMPLETE ~\n')


    # generate lists of tree files to do locus-taxon occupancy guidance_filtering
    with open('supercontig_tre_list.txt', 'w+') as outf:
        for i in to_cat['supercontig']:
            outf.write('{}\n'.format(i))

    with open('exon_tre_list.txt', 'w+') as outf:
        for i in to_cat['inverse']:
            outf.write('{}\n'.format(i))
        for i in to_cat['exon']:
            outf.write('{}\n'.format(i))

    with open('supercontig-exon_tre_list.txt', 'w+') as outf:
        for i in to_cat['supercontig']:
            outf.write('{}\n'.format(i))
        for i in to_cat['inverse']:
            outf.write('{}\n'.format(i))


    # generate lists of alignment files to do locus-taxon occupancy guidance_filtering
    with open('supercontig_aln_list.txt', 'w+') as outf:
        for i in to_aln['supercontig']:
            outf.write('{}\n'.format(i))

    with open('exon_aln_list.txt', 'w+') as outf:
        for i in to_aln['inverse']:
            outf.write('{}\n'.format(i))
        for i in to_aln['exon']:
            outf.write('{}\n'.format(i))

    with open('supercontig-exon_aln_list.txt', 'w+') as outf:
        for i in to_aln['supercontig']:
            outf.write('{}\n'.format(i))
        for i in to_aln['inverse']:
            outf.write('{}\n'.format(i))


    # generate lists of trees
    #with open('{}-{}_supercontigs.txt'.format(loci_filter, loclen_filter), 'w+') as outf:
    #    for f in to_cat[1]:
    #        outf.write('{}\n'.format(f))

    #with open('{}-{}_supercontigs-exons.txt'.format(loci_filter, loclen_filter), 'w+') as outf:
    #    for f in to_cat[1]:
    #        outf.write('{}\n'.format(f))
    #    for f in to_cat[2]:
    #        outf.write('{}\n'.format(f))



    print('\n~ GENERATE AMAS SUMMARIES ~')
    # generate alignments for concatenated analysis
    try:
        os.system('AMAS.py concat -i {} -f fasta -d dna -t supercontig.gene.aln -p supercontig_partitions.txt'.format(' '.join(to_aln['supercontig'])))
        os.system('AMAS.py concat -i {} {} -f fasta -d dna -t supercontig-exon.gene.aln -p supercontig-exon_partitions.txt'.format(' '.join(to_aln['supercontig']), ' '.join(to_aln['inverse'])))
        os.system('AMAS.py concat -i {} -f fasta -d dna -t exon.gene.aln -p exon_partitions.txt'.format(' '.join(to_aln['exon'])))
    except:
        print('did you "conda activate amas"?')
        sys.exit()
    print('~ GENERATE AMAS SUMMARIES COMPLETE ~\n')


    # STATISTICS

    # calculate alignment stats
    # gets alignment length and proportion of parsimony informative characters
    os.system('AMAS.py summary -i *.aln -f fasta -d dna -o summary.txt')

    amas = pandas.read_csv('summary.txt', sep='\t', header=0, index_col=0)
    amas_data = {}
    for index, row in amas.iterrows():
        aln_length = row['Alignment_length']
        pis_perc = row['Proportion_parsimony_informative']
        num_taxa = row['No_of_taxa']

        analysis = index.split('.')[0]
        amas_data[analysis] = [num_taxa, aln_length, pis_perc]

        taxa = []
        with open(index, 'r') as inf:
            for record in SeqIO.parse(inf, 'fasta'):
                taxa.append(record.id)

        taxa = list(set(taxa))
        amas_data[analysis].append(taxa)



    # SUPERCONTIGS
    print('\n\n~~ SUPERCONTIGS')
    print('\n{} total loci in supercontig analyses: {}'.format(len(loci['supercontig']), ','.join(loci['supercontig'])))
    print('\n{} taxa included in supercontig analysis: ({}) {}\n'.format(amas_data['supercontig'][0], len(amas_data['supercontig'][3]), ','.join(amas_data['supercontig'][3])))
    print('\nin the supercontig-only concatenated analysis:')
    print('\talignment length = {}'.format(amas_data['supercontig'][1]))
    print('\tproportion PIS = {}'.format(amas_data['supercontig'][2]))

    # SUPERCONTIGS+EXONS
    print('\n\n~~ SUPERCONTIGS+EXONS')
    print('\n{} exons added: {}'.format(len(loci['inverse']), ','.join(loci['inverse'])))
    print('\n{} total loci in supercontig+exon analyses: {},{}'.format(len(loci['supercontig']) + len(loci['inverse']), ','.join(loci['supercontig']), ','.join(loci['inverse'])))
    print('\n{} taxa included in supercontig+exon analyses: ({}) {}'.format(amas_data['supercontig-exon'][0], len(amas_data['supercontig-exon'][3]), ','.join(amas_data['supercontig-exon'][3])))
    print('\nin the supercontig-exon concatenated analysis:')
    print('\talignment length = {}'.format(amas_data['supercontig-exon'][1]))
    print('\tproportion PIS = {}'.format(amas_data['supercontig-exon'][2]))

    # EXONS
    print('\n\n~~ EXONS')
    print('\n{} total loci in exon analysis: {}'.format(len(loci['exon']), ','.join(loci['exon'])))
    print('\n{} taxa included in exon analyses: ({}) {}'.format(amas_data['exon'][0], len(amas_data['exon'][3]), ','.join(amas_data['exon'][3])))
    print('\nin the exon-only concatenated analysis:')
    print('\talignment length = {}'.format(amas_data['exon'][1]))
    print('\tproportion PIS = {}'.format(amas_data['exon'][2]))

    print('\ntotal loci to avoid: {}'.format(len(avoid)))
    print('{} loci flagged (supercontig): {}'.format(len(flagged_loci), ','.join(flagged_loci)))
    print('{} loci ignored (exon): {}'.format(len(ignored_loci), ','.join(ignored_loci)))


    # modify partitions files
    os.system("sed -i 's/p/DNA,p/' *partitions*")



if __name__ == '__main__':
    main()
