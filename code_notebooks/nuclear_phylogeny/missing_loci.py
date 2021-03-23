#!/usr/bin/python

import sys
import os
import argparse
import pandas
from Bio import SeqIO
from amas import AMAS

def concatenate_alignments(to_aln, loci, flagged_loci, ignored_loci, avoid):

    print('\n~ CONCATENATE ALIGNMENTS ~')
    
    s_aln = AMAS.MetaAlignment(in_files=to_aln['supercontig'], data_type='dna', in_format='fasta', cores=1)
    se_files = to_aln['supercontig'] + to_aln['inverse']
    se_aln = AMAS.MetaAlignment(in_files=se_files, data_type='dna', in_format='fasta', cores=1)
    e_aln = AMAS.MetaAlignment(in_files=to_aln['exon'], data_type='dna', in_format='fasta', cores=1)

    s_parsed = s_aln.get_parsed_alignments()
    s_tuple = s_aln.get_concatenated(s_parsed)
    s_concat = s_tuple[0]
    s_partitions = s_tuple[1]
    with open('supercontig.gene.aln', 'w+') as outf:
        outf.write(s_aln.print_fasta(s_concat))
    with open('supercontig_partitions.txt', 'w+') as outf:
        outf.write(s_aln.print_raxml_partitions('dna'))

    se_parsed = se_aln.get_parsed_alignments()
    se_tuple = se_aln.get_concatenated(se_parsed)
    se_concat = se_tuple[0]
    se_partitions = se_tuple[1]
    with open('supercontig-exon.gene.aln', 'w+') as outf:
        outf.write(se_aln.print_fasta(se_concat))
    with open('supercontig-exon_partitions.txt', 'w+') as outf:
        outf.write(se_aln.print_raxml_partitions('dna'))

    e_parsed = e_aln.get_parsed_alignments()
    e_tuple = e_aln.get_concatenated(e_parsed)
    e_concat = e_tuple[0]
    e_partitions = e_tuple[1]
    with open('exon.gene.aln', 'w+') as outf:
        outf.write(e_aln.print_fasta(e_concat))
    with open('exon_partitions.txt', 'w+') as outf:
        outf.write(e_aln.print_raxml_partitions('dna'))

    print('~ CONCATENATE ALIGNMENTS COMPLETE ~\n')



    # STATISTICS

    print('\n~ GENERATE ALIGNMENT STATISTICS ~')
    # calculate alignment stats
    # gets alignment length and proportion of parsimony informative characters
    alns = []
    for i in os.listdir(os.getcwd()):
        if i.endswith('aln'):
            alns.append(i)
    summary_alns = AMAS.MetaAlignment(in_files=alns, data_type='dna', in_format='fasta', cores=1)
    summary_alns.write_summaries('summary.txt')

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



def generate_aln_list(to_aln):
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
            
            
            
def generate_tre_list(to_cat):
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



def concatenate_gene_trees(to_cat):
    print('\n~ CONCATENATE GENE TREES ~')
    os.system('cat {} > supercontig.gene.tre'.format(' '.join(to_cat['supercontig'])))
    os.system('cat {} {} > supercontig-exon.gene.tre'.format(' '.join(to_cat['supercontig']), ' '.join(to_cat['inverse'])))
    os.system('cat {} > exon.gene.tre'.format(' '.join(to_cat['exon'])))
    print('~ CONCATENATE GENE TREES COMPLETE ~\n')



def only_aligments(avoid, supercontig_aln, exon_aln, bootstrap):
    
    final_loci = []
    flagged_loci = []
    ignored_loci = []
    
    for f in os.listdir(supercontig_aln):
        if f.endswith('fasta'):
            locus = f.split('.')[0]
            if locus not in final_loci:
                final_loci.append(locus)
    
    for f in os.listdir(exon_aln):
        if f.endswith('fasta'):
            locus = f.split('.')[0]
            if locus not in final_loci:
                final_loci.append(locus)
                
    print(len(final_loci))
            
    # first take all supercontigs (supercontig-only)
    # second take all exons for which there were no supercontig dataset (supercontig+exon)
    # third take all exons (exon-only)
    
    to_aln = {'supercontig':[], 'inverse':[], 'exon':[]} # fasta files to concatenate
    loci = {'supercontig':[], 'inverse':[], 'exon':[]} # which loci are being included?
    
    for f in os.listdir(supercontig_aln):
        if f.endswith('.fasta'):
            locus = f.split('.')[0]

            if locus not in avoid:
                if locus in final_loci:
                    to_aln['supercontig'].append(os.path.join(supercontig_aln, f))
                    loci['supercontig'].append(locus)
            else:
                flagged_loci.append(locus)
                    
    
    for f in os.listdir(exon_aln):
        if f.endswith('.fasta'):
            locus = f.split('.')[0]

            if locus not in avoid:
                if locus in final_loci:
                    if locus not in loci['supercontig']:# exon locus not yet accounted for by supercontig data
                        to_aln['inverse'].append(os.path.join(exon_aln, f))
                        loci['inverse'].append(locus)
                    to_aln['exon'].append(os.path.join(exon_aln, f))
                    loci['exon'].append(locus)
            else:
                ignored_loci.append(locus)

                
    # generate lists of alignment files to do locus-taxon occupancy guidance_filtering
    generate_aln_list(to_aln)

    # concatenate alignements and generate amas summaries
    concatenate_alignments(to_aln, loci, flagged_loci, ignored_loci, avoid)


    
    
def only_trees(avoid, supercontig_tre, exon_tre, bootstrap):
    
    
    final_loci = []
    flagged_loci = []
    ignored_loci = []
    
    for f in os.listdir(supercontig_tre):
        if bootstrap:
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
        if bootstrap:
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
    loci = {'supercontig':[], 'inverse':[], 'exon':[]} # which loci are being included?
    # NOTE: taxa are not being accounted for yet -- should be the same between coalescent/concat analyses, so we'll get those later

    for f in os.listdir(supercontig_tre):
        if bootstrap:
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
                    
    for f in os.listdir(exon_tre):
        if bootstrap:
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
                    
                    
    # generate files of trees for coalescent analysis
    concatenate_gene_trees(to_cat)

    # generate lists of tree files to do locus-taxon occupancy guidance_filtering
    generate_tre_list(to_cat)

    print("\n\nyou have provided only trees. unfortunately, I cannot generate alignment statistics from these.")
    print("missing_loci.py complete.\n\n")



    
    
    
def missing_loci(avoid, supercontig_aln, exon_aln, supercontig_tre, exon_tre, bootstrap):
    # get which loci are inverse between the supercontig/exon trees and alignments
    # (some loci will have dropped out due to RAxML only taking loci with >=4 taxa)
    # ensures consistency between coalescent and concatenated analyses
    
    final_loci = []
    flagged_loci = []
    ignored_loci = []
    
    for f in os.listdir(supercontig_tre):
        if bootstrap:
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
        if bootstrap:
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
        if bootstrap:
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
            locus = f.split('.')[0][:4]

            if locus in loci['supercontig']:
                if locus not in avoid:
                    if locus in final_loci:
                        to_aln['supercontig'].append(os.path.join(supercontig_aln, f))

   
    for f in os.listdir(exon_tre):
        if bootstrap:
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
            locus = f.split('.')[0][:4]

            if locus in loci['inverse']:
                to_aln['inverse'].append(os.path.join(exon_aln, f))
            
            if locus in loci['exon']:
                to_aln['exon'].append(os.path.join(exon_aln, f))

    print(to_cat['exon'], to_aln['exon'])
    # generate files of trees for coalescent analysis
    concatenate_gene_trees(to_cat)

    # generate lists of tree files to do locus-taxon occupancy guidance_filtering
    generate_tre_list(to_cat)

    # generate lists of alignment files to do locus-taxon occupancy guidance_filtering
    generate_aln_list(to_aln)

    # concatenate alignements and generate amas summaries
    #concatenate_alignments(to_aln, loci, flagged_loci, ignored_loci, avoid)

    
    


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
    
    if args.supercontig_aln:
        supercontig_aln = os.path.abspath(args.supercontig_aln)
    if args.supercontig_tre:
        supercontig_tre = os.path.abspath(args.supercontig_tre)
    if args.exon_aln:
        exon_aln = os.path.abspath(args.exon_aln)
    if args.exon_tre:
        exon_tre = os.path.abspath(args.exon_tre)
        
    bootstrap = args.bootstrap


    # these are loci for which the recovered exons did reach a certain percentage of the target locus
    # if found in supercontig dataset, flagged but included due to flanking regions providing additional information
    # if found in exon dataset, ignored.
    avoid = []
    if args.avoid_loci:
        with open(args.avoid_loci, 'r') as inf:
            for line in inf:
                avoid.append(line.strip())
    
    
    # combine loci
    # only alignments
    if (args.supercontig_aln and args.exon_aln) and not (args.supercontig_tre and args.exon_tre):
        only_aligments(avoid, supercontig_aln, exon_aln, bootstrap)
    
    elif (args.supercontig_tre and args.exon_tre) and not (args.supercontig_aln and args.exon_aln):
        only_trees(avoid, supercontig_tre, exon_tre, bootstrap)
        
    else:
        missing_loci(avoid, supercontig_aln, exon_aln, supercontig_tre, exon_tre, bootstrap)



    



if __name__ == '__main__':
    main()
