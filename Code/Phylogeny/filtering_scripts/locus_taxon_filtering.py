#!/usr/bin/python

# 06-May-2020, Aaron K. Lee

# filter out loci on the basis of how many taxa are represented

# from the input list: identifies the tree / alignment with the largest number of tips / taxa and sets that as MAX
# generates concatnated .aln / .tre files with all input alignments / trees that meet a threshold of ( MAX * 0.xx ) taxa
# xx = 25, 50, 65, 80


# input: a text file list of *trees* or alignments
# output: four .tre or .aln files


import sys, os
from Bio import SeqIO
import pandas
import argparse
import dendropy


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--list_file', help='a list of alignment or tree files, to be parsed by the script, do not use this if using -d *')
    parser.add_argument('-t', '--list_type', help='"aln" or "tre" to indicate the kinds of files in the list file, do not use this if using -d *')
    parser.add_argument('-d', '--collapsed_tree_directory', help='path to directory containing collapsed trees, do not use this if using -l * -t *')
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    if args.list_file:
        listf = os.path.abspath(args.list_file)
        listtype = args.list_type
        if listtype == 'tre':
            filter_tre(listf)
        elif listtype == 'aln':
            filter_aln(listf)

    elif args.collapsed_tree_directory:
        treedir = os.path.abspath(args.collapsed_tree_directory)
        filter_collapsed(treedir)




def filter_tre(listf):
    tre_files = []

    with open(listf, 'r') as inf:
        for line in inf:
            tre_files.append(line.strip())

    locus_taxa = {}
    max = 0
    for f in tre_files:
        with open(f, 'r') as inf:
            locus = f.split('.')[1]
            num = inf.readlines()[0].count(',') + 1
            locus_taxa[locus] = num
            if num > max:
                max = num


    t10, t25, t50, t65, t80 = [], [], [], [], []
    rejected = []
    for l, n in locus_taxa.items():
        if n > int(0.80*max):
            t10.append(l)
            t25.append(l)
            t50.append(l)
            t65.append(l)
            t80.append(l)
        elif n > int(0.65*max):
            t10.append(l)
            t25.append(l)
            t50.append(l)
            t65.append(l)
        elif n > int(0.5*max):
            t10.append(l)
            t25.append(l)
            t50.append(l)
        elif n > int(0.25*max):
            t10.append(l)
            t25.append(l)
        elif n > int(0.10*max):
            t10.append(l)
        else:
            rejected.append(l)

    t25_f = '{}-{}_loc_tax_num.txt'.format(listf.split('.')[0], int(0.25*max))
    t50_f = '{}-{}_loc_tax_num.txt'.format(listf.split('.')[0], int(0.50*max))
    t65_f = '{}-{}_loc_tax_num.txt'.format(listf.split('.')[0], int(0.65*max))
    t80_f = '{}-{}_loc_tax_num.txt'.format(listf.split('.')[0], int(0.80*max))

    t10_t = '{}-{}_loc_tax_num.gene.tre'.format(listf.split('/')[-1].split('.')[0].split('_')[0], '10max')
    t25_t = '{}-{}_loc_tax_num.gene.tre'.format(listf.split('/')[-1].split('.')[0].split('_')[0], '25max')
    t50_t = '{}-{}_loc_tax_num.gene.tre'.format(listf.split('/')[-1].split('.')[0].split('_')[0], '50max')
    t65_t = '{}-{}_loc_tax_num.gene.tre'.format(listf.split('/')[-1].split('.')[0].split('_')[0], '65max')
    t80_t = '{}-{}_loc_tax_num.gene.tre'.format(listf.split('/')[-1].split('.')[0].split('_')[0], '80max')


    t10_trees = []
    for f in tre_files:
        locus = f.split('.')[1]
        if locus in t10:
            t10_trees.append(f)
    os.system('cat {} > {}'.format(' '.join(t10_trees), t10_t))


    t25_trees = []
    for f in tre_files:
        locus = f.split('.')[1]
        if locus in t25:
            t25_trees.append(f)
    os.system('cat {} > {}'.format(' '.join(t25_trees), t25_t))

    t50_trees = []
    for f in tre_files:
        locus = f.split('.')[1]
        if locus in t50:
            t50_trees.append(f)
    os.system('cat {} > {}'.format(' '.join(t50_trees), t50_t))

    t65_trees = []
    for f in tre_files:
        locus = f.split('.')[1]
        if locus in t65:
            t65_trees.append(f)
    os.system('cat {} > {}'.format(' '.join(t65_trees), t65_t))

    t80_trees = []
    for f in tre_files:
        locus = f.split('.')[1]
        if locus in t80:
            t80_trees.append(f)
    os.system('cat {} > {}'.format(' '.join(t80_trees), t80_t))


    print('started out with {} loci'.format(len(locus_taxa)))

    print('\n\n10% of {} (max num taxa) = at least {} taxa per locus'.format(max, int(0.1*max)))
    print('resulting in {} loci\n'.format(len(t10_trees)))

    print('25% of {} (max num taxa) = at least {} taxa per locus'.format(max, int(0.25*max)))
    print('resulting in {} loci\n'.format(len(t25_trees)))

    print('50% of {} (max num taxa) = at least {} taxa per locus'.format(max, int(0.5*max)))
    print('resulting in {} loci\n'.format(len(t50_trees)))

    print('65% of {} (max num taxa) = at least {} taxa per locus'.format(max, int(0.65*max)))
    print('resulting in {} loci\n'.format(len(t65_trees)))

    print('80% of {} (max num taxa) = at least {} taxa per locus'.format(max, int(0.80*max)))
    print('resulting in {} loci\n\n'.format(len(t80_trees)))

    print('{} loci did not meet percentage thresholds: {}\n\n'.format(len(rejected), ','.join(rejected)))




def filter_aln(listf):
    aln_files = []

    with open(listf, 'r') as inf:
        for line in inf:
            aln_files.append(line.strip())

    locus_taxa = {}
    max = 0
    for f in aln_files:
        with open(f, 'r') as inf:
            locus = f.split('/')[-1][0:4]
            num = 0
            for record in SeqIO.parse(inf, 'fasta'):
                num += 1
            locus_taxa[locus] = num
            if num > max:
                max = num

    a10, a25, a50, a65, a80 = [], [], [], [], []
    rejected = []
    for l, n in locus_taxa.items():
        if n > int(0.80*max):
            a10.append(l)
            a25.append(l)
            a50.append(l)
            a65.append(l)
            a80.append(l)
        elif n > int(0.65*max):
            a10.append(l)
            a25.append(l)
            a50.append(l)
            a65.append(l)
        elif n > int(0.5*max):
            a10.append(l)
            a25.append(l)
            a50.append(l)
        elif n > int(0.25*max):
            a10.append(l)
            a25.append(l)
        elif n > int(0.10*max):
            a10.append(l)
        else:
            rejected.append(l)


    a10_a = '{}-{}_loc_tax_num.aln'.format(listf.split('/')[-1].split('.')[0].split('_')[0], '10max')
    a25_a = '{}-{}_loc_tax_num.aln'.format(listf.split('/')[-1].split('.')[0].split('_')[0], '25max')
    a50_a = '{}-{}_loc_tax_num.aln'.format(listf.split('/')[-1].split('.')[0].split('_')[0], '50max')
    a65_a = '{}-{}_loc_tax_num.aln'.format(listf.split('/')[-1].split('.')[0].split('_')[0], '65max')
    a80_a = '{}-{}_loc_tax_num.aln'.format(listf.split('/')[-1].split('.')[0].split('_')[0], '80max')

    a10_p = '{}-{}_loc_tax_num.partitions.txt'.format(listf.split('/')[-1].split('.')[0].split('_')[0], '10max')
    a25_p = '{}-{}_loc_tax_num.partitions.txt'.format(listf.split('/')[-1].split('.')[0].split('_')[0], '25max')
    a50_p = '{}-{}_loc_tax_num.partitions.txt'.format(listf.split('/')[-1].split('.')[0].split('_')[0], '50max')
    a65_p = '{}-{}_loc_tax_num.partitions.txt'.format(listf.split('/')[-1].split('.')[0].split('_')[0], '65max')
    a80_p = '{}-{}_loc_tax_num.partitions.txt'.format(listf.split('/')[-1].split('.')[0].split('_')[0], '80max')

    a10_aln = []
    for f in aln_files:
        locus = f.split('/')[-1][0:4]
        if locus in a10:
            a10_aln.append(f)
    os.system('AMAS.py concat -i {} -f fasta -d dna -t {} -p {}'.format(' '.join(a10_aln), a10_a, a10_p))

    a25_aln = []
    for f in aln_files:
        locus = f.split('/')[-1][0:4]
        if locus in a25:
            a25_aln.append(f)
    os.system('AMAS.py concat -i {} -f fasta -d dna -t {} -p {}'.format(' '.join(a25_aln), a25_a, a25_p))

    a50_aln = []
    for f in aln_files:
        locus = f.split('/')[-1][0:4]
        if locus in a50:
            a50_aln.append(f)
    os.system('AMAS.py concat -i {} -f fasta -d dna -t {} -p {}'.format(' '.join(a50_aln), a50_a, a50_p))

    a65_aln = []
    for f in aln_files:
        locus = f.split('/')[-1][0:4]
        if locus in a65:
            a65_aln.append(f)
    os.system('AMAS.py concat -i {} -f fasta -d dna -t {} -p {}'.format(' '.join(a65_aln), a65_a, a65_p))

    a80_aln = []
    for f in aln_files:
        locus = f.split('/')[-1][0:4]
        if locus in a80:
            a80_aln.append(f)
    os.system('AMAS.py concat -i {} -f fasta -d dna -t {} -p {}'.format(' '.join(a80_aln), a80_a, a80_p))

    os.system("sed -i 's/p/DNA,p/' *partitions*")


    # calculate alignment statistics
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


    # output results of filtering
    print('started out with {} loci'.format(len(locus_taxa)))

    print('\n\n10% of {} (max num taxa) = at least {} taxa per locus'.format(max, int(0.1*max)))
    print('resulting in {} loci\n'.format(len(a10_aln)))

    print('25% of {} (max num taxa) = at least {} taxa per locus'.format(max, int(0.25*max)))
    print('resulting in {} loci\n'.format(len(a25_aln)))

    print('50% of {} (max num taxa) = at least {} taxa per locus'.format(max, int(0.5*max)))
    print('resulting in {} loci\n'.format(len(a50_aln)))

    print('65% of {} (max num taxa) = at least {} taxa per locus'.format(max, int(0.65*max)))
    print('resulting in {} loci\n'.format(len(a65_aln)))

    print('80% of {} (max num taxa) = at least {} taxa per locus'.format(max, int(0.80*max)))
    print('resulting in {} loci\n\n'.format(len(a80_aln)))

    print('alignment statistics:')
    for i, j in amas_data.items():
        print('{}'.format(i))
        print('\tnum taxa = {}'.format(j[0]))
        print('\talignment length = {}'.format(j[1]))
        print('\tpis_perc = {}\n\n'.format(j[2]))

    print('{} loci did not meet percentage thresholds: {}\n\n'.format(len(rejected), ','.join(rejected)))




def filter_collapsed(treedir):
    data = ['exon', 'supercontig', 'supercontig-exon']
    tree_dict = {i:[] for i in data}

    for file in os.listdir(treedir):
        d = file.split('.')[0]
        if d in data:
            tree_dict[d].append(os.path.join(treedir, file))

    for i, j in tree_dict.items():
        locus_taxa = {}
        max = 0
        for k in j:
            locus = k.split('.')[1].split('_')[0]
            tre = dendropy.Tree.get(path=k, schema='newick')
            num = 0
            for node in tre.seed_node.leaf_iter():
                num += 1

            locus_taxa[locus] = num

            if num > max:
                max = num


        t10, t25, t50, t65, t80 = [], [], [], [], []
        rejected = []
        for l, n in locus_taxa.items():
            if n > int(0.80*max):
                t10.append(l)
                t25.append(l)
                t50.append(l)
                t65.append(l)
                t80.append(l)
            elif n > int(0.65*max):
                t10.append(l)
                t25.append(l)
                t50.append(l)
                t65.append(l)
            elif n > int(0.5*max):
                t10.append(l)
                t25.append(l)
                t50.append(l)
            elif n > int(0.25*max):
                t10.append(l)
                t25.append(l)
            elif n > int(0.10*max):
                t10.append(l)
            else:
                rejected.append(l)


        t10_t = '{}-{}_loc_tax_num.gene.tre'.format(i, '10max')
        t25_t = '{}-{}_loc_tax_num.gene.tre'.format(i, '25max')
        t50_t = '{}-{}_loc_tax_num.gene.tre'.format(i, '50max')
        t65_t = '{}-{}_loc_tax_num.gene.tre'.format(i, '65max')
        t80_t = '{}-{}_loc_tax_num.gene.tre'.format(i, '80max')


        t10_trees = []
        for f in j:
            locus = f.split('.')[1].split('_')[0]
            if locus in t10:
                t10_trees.append(f)
        os.system('cat {} > {}'.format(' '.join(t10_trees), t10_t))


        t25_trees = []
        for f in j:
            locus = f.split('.')[1].split('_')[0]
            if locus in t25:
                t25_trees.append(f)
        os.system('cat {} > {}'.format(' '.join(t25_trees), t25_t))

        t50_trees = []
        for f in j:
            locus = f.split('.')[1].split('_')[0]
            if locus in t50:
                t50_trees.append(f)
        os.system('cat {} > {}'.format(' '.join(t50_trees), t50_t))

        t65_trees = []
        for f in j:
            locus = f.split('.')[1].split('_')[0]
            if locus in t65:
                t65_trees.append(f)
        os.system('cat {} > {}'.format(' '.join(t65_trees), t65_t))

        t80_trees = []
        for f in j:
            locus = f.split('.')[1].split('_')[0]
            if locus in t80:
                t80_trees.append(f)
        os.system('cat {} > {}'.format(' '.join(t80_trees), t80_t))


        print(i)
        print('started out with {} loci'.format(len(locus_taxa)))

        print('\n10% of {} (max num taxa) = at least {} taxa per locus'.format(max, int(0.1*max)))
        print('resulting in {} loci\n'.format(len(t10_trees)))

        print('25% of {} (max num taxa) = at least {} taxa per locus'.format(max, int(0.25*max)))
        print('resulting in {} loci\n'.format(len(t25_trees)))

        print('50% of {} (max num taxa) = at least {} taxa per locus'.format(max, int(0.5*max)))
        print('resulting in {} loci\n'.format(len(t50_trees)))

        print('65% of {} (max num taxa) = at least {} taxa per locus'.format(max, int(0.65*max)))
        print('resulting in {} loci\n'.format(len(t65_trees)))

        print('80% of {} (max num taxa) = at least {} taxa per locus'.format(max, int(0.80*max)))
        print('resulting in {} loci\n\n'.format(len(t80_trees)))

        print('{} loci did not meet percentage thresholds: {}\n\n'.format(len(rejected), ','.join(rejected)))




if __name__ == '__main__':
    main()
