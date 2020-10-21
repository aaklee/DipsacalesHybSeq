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
            region = ''
            for index, row in ann_subunits_df.iterrows():
                if 'junction' in ann_df['Name']:
                    continue
                else:
                    name = row['Name'].lower()
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


def parse_intergenic_regions(annot, refspecies, subunits):
    # parse just one species annotation into intergenic regions
    print(refspecies)
    intergenic = {}

    currmin = 0
    currmax = 0
    curr_intergenic = []

    for i in range(len(annot)):
        j = i + 1

        currgene = annot[i]
        if '>' in currgene[2]:
            continue
        currname = currgene[0]
        if currmin == 0 and currmax == 0:
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

        if currmax > nextmin:
            currmax = nextmax
            curr_intergenic.append(nextname)
            continue

        itmin = currmax
        itmax = nextmin
        curr_intergenic.extend([currname, nextname])
        new_region = '_'.join(curr_intergenic)


        # check if intergenic region overlaps with subunit boundaries
        overlap = False
        for x in range(len(subunits[refspecies])):
            xlower = int(subunits[refspecies][x][1])
            xupper = int(subunits[refspecies][x][2])
            xregion = subunits[refspecies][x][0]

            # the intergenic region is safely nestled in a cp subunit
            if itmin > xlower and itmax < xupper:
                print('{} in {}'.format(new_region, xregion))
                break


            # the intergenic region is at the junction of two subunits
            y = x + 1
            # at the circular part of the plastid
            if y == len(subunits[refspecies]):
                ylower = int(subunits[refspecies][0][1])
                yupper = int(subunits[refspecies][0][2])
                yregion = subunits[refspecies][0][0]

            # at any other junction
            else:
                ylower = int(subunits[refspecies][y][1])
                yupper = int(subunits[refspecies][y][2])
                yregion = subunits[refspecies][y][0]

            if itmin > xlower and itmax > ylower:
                # doesn't belong to x subunit
                if itmin < xupper and itmax < yupper:
                    # at a junction
                    print('WARNING: {} in between {} and {}'.format(new_region, xregion, yregion))
                    overlap = True


        if not overlap:
            intergenic[new_region] = (itmin, itmax)
            #print(new_region, intergenic[new_region])


        # reset
        currmin, currmax = 0, 0
        curr_intergenic = []


    print()
    return intergenic



def parse_fasta(refspecies, intergenic, fastas):
    seqs = fastas[refspecies]
    for s, t in seqs.items():
        final_seqs = {}
        seq_array = list(t)

        for i, j in intergenic.items():
            new_intergenic = '.'.join([s, i]) # species.region
            start = int(j[0]-1)
            end = int(j[1])

            currseq = ''.join(seq_array[start:end])
            # clean out ambiguous bp and gaps
            currseq = currseq.replace('?', '')
            currseq = currseq.replace('-', '')

            if len(currseq) > 0:
                final_seqs[new_intergenic] = currseq

        newf = '{}_intergenic_spacers.fasta'.format(s)
        with open(newf, 'w+') as outf:
            for x, y in final_seqs.items():
                outf.write('>{}\n'.format(x))
                outf.write('{}\n'.format(y))




def extract_intergenic_regions(aln, ann):
    fastas = parse_aln(aln)
    convert_csv(ann)
    annotations, subunits = parse_ann(ann)

    # parse out intergenic regions
    for refspecies, annot in annotations.items():
        intergenic = parse_intergenic_regions(annot, refspecies, subunits)
        parse_fasta(refspecies, intergenic, fastas)




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

    extract_intergenic_regions(aln, ann)

if __name__ == '__main__':
    main()
