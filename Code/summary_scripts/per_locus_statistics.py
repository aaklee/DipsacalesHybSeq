#!/usr/bin/python

# aaron lee, 15-june-2020
# generate statistics and plots given a list of trees and a list of alignments
# currently:
# - %PIS histogram per locus
# - %resolution (% nodes supported by at least 50% of gene trees)

import sys
import os
from matplotlib import pyplot
import argparse
import pandas
import numpy
import dendropy

OUTGROUPS = ['Ilex crenata', 'Paracryphia alticola']

def reroot(tree):
    root_node = dendropy.Node()

    for i in OUTGROUPS:
        print(i)
        root_node = tree.find_node_with_taxon_label(i)
        if root_node:
            break

    print(root_node)
    if root_node:
        tree.reroot_at_edge(root_node.edge)
    return tree


def calculate_supporting_bipartitions(tree_list, species_tree):
    gene_trees = []
    genes = []
    with open(tree_list, 'r') as inf:
        for line in inf:
            gene_trees.append(line.strip())
            gene = line.strip().split('/')[-1].split('.')[1]
            genes.append(gene)

    trees = dendropy.TreeList()
    # trees[0] is always the species tree
    trees.append(dendropy.Tree.get(path=species_tree, schema='newick')) #, rooting='force-rooted'))
    #trees[0] = reroot(trees[0])

    for i in range(len(gene_trees)):
        trees.append(dendropy.Tree.get(path=gene_trees[i], schema='newick')) #, rooting='force-rooted'))
        #trees[i+1] = reroot(trees[i+1])

    sp_splits = trees[0].encode_bipartitions()
    fin_topology = {x:[] for x in [y.leafset_as_bitstring() for y in sp_splits]}

    alt_bipartitions = {x:[] for x in genes}


    for i in range(1, len(trees)):
        for j in trees[i].encode_bipartitions():
            k_bipts = []
            for k in sp_splits:
                if j.is_leafset_nested_within(k):
                    fin_topology[k.leafset_as_bitstring()].append(genes[i-1])

    # summarize
    # bitwise bipartitions; how many/what percent of gene trees have a bipartition that is nested within?
    for k, v in fin_topology.items():
        if str(k).count('1') != 1:
            vs = set(v)
            print('{}\t{}/{}\t{}'.format(k, len(vs), len(genes), float(len(vs)/len(genes))))




def calculate_perc_resolution(tree_list, species_tree):
    support_threshold = 50

    gene_trees = []
    genes = []
    with open(tree_list, 'r') as inf:
        for line in inf:
            gene_trees.append(line.strip())
            gene = line.strip().split('/')[-1].split('.')[1]
            genes.append(gene)

    trees = dendropy.TreeList()
    # trees[0] is always the species tree
    trees.append(dendropy.Tree.get(path=species_tree, schema='newick'))

    for i in range(len(gene_trees)):
        trees.append(dendropy.Tree.get(path=gene_trees[i], schema='newick'))

    # percent resolution for the species tree
    supported_nodes = 0
    possible_nodes = len(trees[0].nodes())
    if 'RAx' in species_tree:
        for edge in trees[0].postorder_edge_iter():
            if edge.head_node.label and int(edge.head_node.label) >= support_threshold:
                supported_nodes += 1
    else:
        for edge in trees[0].postorder_edge_iter():
            if edge.head_node.label and int(float(edge.head_node.label) * 100) >= support_threshold:
                supported_nodes += 1

    print('\nspecies tree % resolution:')
    print('\t{}/{} = {:.3f}\n'.format(supported_nodes, possible_nodes, supported_nodes/possible_nodes))

    # percent resolution for gene trees

    print('gene tree percent resolution:')
    for i in range(1, len(trees)):
        supported_nodes = 0
        possible_nodes = len(trees[i].nodes())

        for edge in trees[i].postorder_edge_iter():
            if edge.head_node.label and int(edge.head_node.label) >= support_threshold:
                supported_nodes += 1

        print('{}\t{}/{} = {:.3f}'.format(genes[i-1], supported_nodes, possible_nodes, supported_nodes/possible_nodes))




def generate_pis_plot(aln_list):
    aln = []
    with open(aln_list, 'r') as inf:
        for line in inf:
            aln.append(line.strip())

    #name_parts = aln_list.split('/')
    #name_parts[-3] = int(float(name_parts[-3]) * 100)

    # run AMAS summary to get stats
    #sum_fn = 'locus_summary.{}.{}.{}.txt'.format(name_parts[-4], name_parts[-3], name_parts[-2])
    #print(sum_fn)
    #os.system('AMAS.py summary -i {} -f fasta -d dna -o {}'.format(' '.join(aln), sum_fn))

    os.system('AMAS.py summary -i {} -f fasta -d dna -o locus_summary.txt'.format(' '.join(aln)))

    # convert stats to usable format as pandas dataframe
    #amas = pandas.read_csv(sum_fn, sep='\t', header=0, index_col=0)

    amas = pandas.read_csv('locus_summary.txt', sep='\t', header=0, index_col=0)
    # generate bar chart
    amas = amas.sort_values(by=['Proportion_parsimony_informative'])
    x = amas.index.to_list()
    y = amas['Proportion_parsimony_informative'].to_list()
    b = amas['Alignment_length'].to_list()
    j = amas['No_of_taxa'].to_list()

    #fn = 'pis_barplot.{}.{}.{}.pdf'.format(name_parts[-4], name_parts[-3], name_parts[-2])
    fig, (axy, axb, axj) = pyplot.subplots(3)
    axy.bar(x, y, align='center')
    axy.set_ylabel('%PIS')
    axy.tick_params(axis='x', labelbottom=False)

    axb.bar(x, b, align='center')
    axb.axhline(y=500, linewidth=0.5, color='r')
    axb.set_ylabel('aln_len')
    axb.set_xlabel('locus')
    axb.tick_params(axis='x', labelbottom=False)

    axj.bar(x, j, align='center')
    axj.set_ylabel('no_taxa')
    axj.tick_params(axis='x', labelbottom=False)

    fig.savefig('pis_barplot.pdf')
    #fig.savefig(fn)


#    fn = 'pis_barplot.{}.{}.{}.pdf'.format(name_parts[-4], name_parts[-3], name_parts[-2])
#    fig, ax = pyplot.subplots()
#    ax.bar(x, y, align='center')
#    ax.set_xlabel('locus')
#    ax.set_ylabel('Proportion_parsimony_informative')
#    ax.tick_params(axis='x', labelbottom=False)
#    fig.savefig(fn)


    # generate another bar chart
#    amas = amas.sort_values(by=['Alignment_length'])
#    x = amas.index.to_list()
#    y = amas['Alignment_length'].to_list()

#    name_parts = aln_list.split('/')
#    fn = 'length_barplot.{}.{}.{}.pdf'.format(name_parts[-4], name_parts[-3], name_parts[-2])
#    fig, ax = pyplot.subplots()
#    ax.bar(x, y, align='center')
#    ax.set_xlabel('locus')
#    ax.set_ylabel('Alignment_length')
#    ax.tick_params(axis='x', labelbottom=False)
#    fig.savefig(fn)

    # generate binned barplot, 20 bins
#    max = amas['Proportion_parsimony_informative'].max()
#    min = amas['Proportion_parsimony_informative'].min()
#    partition = max/20

#    x, y = [], []
#    for i in range(1,21):
#        cutoff1 = (i-1)*partition
#        cutoff2 = i*partition
#        x.append('{:.3f}-<{:.3f}'.format(cutoff1, cutoff2))

#        count = 0
#        for j in amas['Proportion_parsimony_informative']:
#            if j >= cutoff1 and j < cutoff2:
#                count += 1

#        y.append(count)

#    fig, ax = pyplot.subplots()
#    ax.bar(x, y, align='center')
#    ax.set_xlabel('% PIS range')
#    ax.set_ylabel('counts')
#    ax.tick_params(axis='x', labelsize=3, labelrotation=90)
#    fig.savefig('pis_binned_barplot.pdf')





def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tree_list', help='file containing a list of paths to RAxML tree outputs')
    parser.add_argument('-a', '--aln_list', help='file containing a list of paths to alignments in FASTA format')
    parser.add_argument('-g', '--graph_pis', action='store_true', help='include if you want to generate percent PIS histogram per locus')
    parser.add_argument('-r', '--perc_resolution', action='store_true', help='requires "-s", include if you want to calculate gene-wise percent resolution')
    parser.add_argument('-s', '--species_tree', help='path to species tree, for use in "-r"')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()
    if args.perc_resolution:
        if not args.species_tree:
            print('missing "-s" species tree')
            parser.print_help()
            sys.exit()

    # parse arguments
    if args.tree_list:
        tree_list = os.path.abspath(args.tree_list)
    if args.aln_list:
        aln_list = os.path.abspath(args.aln_list)
    graph_pis = args.graph_pis
    perc_resolution = args.perc_resolution

    if perc_resolution:
        species_tree = os.path.abspath(args.species_tree)
        calculate_perc_resolution(tree_list, species_tree)

    if graph_pis:
        generate_pis_plot(aln_list)




if __name__ == '__main__':
    main()
