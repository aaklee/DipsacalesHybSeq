#!/usr/bin/python

import sys
import os
import argparse
from Bio import SeqIO
import pandas

def parse_aln(aln):
    # break down fasta alignment files
    # {species in file name: {reference assembly: seq, ...}}

    fastas = {}
    for i in os.listdir(aln):
        j = i.split('.')[0]
        ref = '_'.join(j.split('_')[:2])
        if ref not in fastas.keys():
            fastas[ref] = {}

        for record in SeqIO.parse(os.path.join(aln, i), 'fasta'):

            if '-' not in record.id:
                species = '_'.join(record.id.split('_')[:2])
                fastas[ref][species] = record.seq


        #print(ref, fastas[ref].keys())

    return fastas


def convert_csv(ann):
    # convert csv to tsv, string to int problem!! very annoying.

    for i in os.listdir(ann):
        if i.endswith('.csv'):
            ann_df = pandas.read_csv(os.path.join(ann, i), sep=',', header=0)
            newf = i.split('.')[0] + '.tsv'
            ann_df.to_csv(os.path.join(ann, newf), sep='\t', header=True, index=False)

    for i in os.listdir(ann):
        if i.endswith('.tsv'):
            os.system('sed -i "s/,//g" {}'.format(os.path.join(ann, i)))

    return


def parse_ann(ann):
    # break down annotations
    # {species in file name: {gene blocks: x-y, ... }}

    annotations = {}
    subunits = {}

    for i in os.listdir(ann):
        if i.endswith('.tsv'):
            species = i.split('.')[0]
            if species not in annotations.keys():
                annotations[species] = []
                subunits[species] = []

            #print(i)
            ann_df = pandas.read_csv(os.path.join(ann, i), sep='\t', header=0)
            ann_df = ann_df.loc[ann_df['Sequence Name'] == 'Consensus']

            # get subunits
            ann_subunits_df = ann_df.loc[(ann_df['Type'] == 'misc_feature') | (ann_df['Type'] == 'repeat_region')]

            # subset only gene regions
            ann_df = ann_df.loc[ann_df['Type'] == 'gene']
            ann_df['Minimum'] = pandas.to_numeric(ann_df['Minimum'])
            ann_df = ann_df.sort_values(by='Minimum')
            #print(ann_df)

            for index, row in ann_df.iterrows():
                gene = row['Name'].split(' ')[0]
                annotations[species].append((gene, row['Minimum'], row['Maximum']))

            # get bases where subunits change over
            #print(ann_subunits_df)
            ann_subunits_df = ann_subunits_df.sort_values(by='Minimum')
            for index, row in ann_subunits_df.iterrows():
                if 'junction' in ann_df['Name']:
                    continue
                else:
                    name = row['Name'].lower()
                    region = ''
                    if 'irb' in name:
                        region = 'IRB'
                    elif 'ira' in name:
                        region = 'IRA'
                    elif 'lsc' in name:
                        region = 'LSC'
                    elif 'ssc' in name:
                        region = 'SSC'
                    subunits[species].append((region, row['Minimum'], row['Maximum']))

            #print(subunits)


    return annotations, subunits


def parse_genic_regions(annot, refspecies, subunits):
    # parse just one species annotation into genic regions
    print(refspecies)
    genes = {}

    gmin = 0
    gmax = 0
    curr_gene = []

    for i in range(len(annot)):
        j = i + 1

        currgene = annot[i]
        if '>' in currgene[2]:
            continue
        currname = currgene[0]
        curr_gene.append(currname)
        currmin = int(currgene[1])
        currmax = int(currgene[2])

        # when the gene is at the end of the linear cp genome
        # compare back to the very first gene
        if j == len(annot):
            nextgene = annot[0]
        else:
            nextgene = annot[j]
        if '>' in nextgene[2]:
            continue
        nextname = nextgene[0]
        nextmin = int(nextgene[1])
        nextmax = int(nextgene[2])

        # if the two gene regions overlap
        if currmax > nextmin:
            gmax = nextmax
            gmin = currmin
            curr_gene.append(nextname)
            continue

        if gmin == 0:
            gmin = currmin
        if gmax == 0:
            gmax = currmax
        print(curr_gene, gmin, gmax)
        new_region = '_'.join(curr_gene) # if the gene region spans >1 gene


        genes[new_region] = (gmin, gmax)
        #print(new_region, genes[new_region])


        # reset
        gmin, gmax = 0, 0
        curr_gene = []


    print()
    return genes



def parse_fasta(refspecies, genes, fastas):
    seqs = fastas[refspecies]
    for s, t in seqs.items():
        final_seqs = {}
        seq_array = list(t)

        print(s, len(seq_array))

        for i, j in genes.items():
            new_genes = '.'.join([s, i]) # species.region
            start = int(j[0]-1)
            end = int(j[1])

            currseq = ''.join(seq_array[start:end])
            # clean out ambiguous bp and gaps
            cleanseq = currseq.replace('?', '')
            cleanseq = cleanseq.replace('-', '')

            if len(cleanseq) > 0:
                print(new_genes, 'added')
                final_seqs[new_genes] = cleanseq
                continue

            print(new_genes, 'NOT added', start, end)
            print(currseq)

        newf = '{}_genic_regions.fasta'.format(s)
        with open(newf, 'w+') as outf:
            for x, y in final_seqs.items():
                outf.write('>{}\n'.format(x))
                outf.write('{}\n'.format(y))





def extract_genic_regions(aln, ann):
    fastas = parse_aln(aln)
    convert_csv(ann)
    annotations, subunits = parse_ann(ann)

    # parse out intergenic regions
    for refspecies, annot in annotations.items():
        genes = parse_genic_regions(annot, refspecies, subunits)
        parse_fasta(refspecies, genes, fastas)




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-aln', help='reference alignments directory, fasta files')
    parser.add_argument('-ann', help='annotations directory, csv files')

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()
    aln = os.path.abspath(args.aln)
    ann = os.path.abspath(args.ann)

    extract_genic_regions(aln, ann)

if __name__ == '__main__':
    main()
