#!/usr/bin/python

# 12-may-2020, aaron lee
# pipeline framework to filter hybpiper results

# input: path to directory containing hybpiper output
# output: everything you ever wanted

import sys
import os
import argparse
from Bio import SeqIO

def retrieve_sequences(sampledir, taxon_pool, STEP, convert_names): # shamelessly copied & modified from hybpiper/retrieve_sequences.py
    retrieve = ['dna', 'intron', 'supercontig', 'aa']
    for r in retrieve:
        if r == 'dna':
            seq_dir = "FNA"
        elif r == 'intron':
            seq_dir = 'intron'
            filename = 'introns'
        elif r == 'supercontig':
            seq_dir = 'intron'
            filename = 'supercontig'
        elif r == 'aa':
            seq_dir = 'FAA'

        #Use gene names parsed from a bait file.
        baitfile = os.path.abspath('/projects/clement-lab/5-HybSeq/2-hybpiper_Dipsacales/Dipsacales_Baits_Project/Angiosperms353_targetSequences.fasta')
        target_genes_dict = SeqIO.to_dict(SeqIO.parse(baitfile,'fasta'))
        target_genes = list(set([x.id.split('-')[-1] for x in SeqIO.parse(baitfile,'fasta')]))

        sample_names = [x for x in os.listdir(sampledir) if (os.path.isdir(os.path.join(sampledir,x)) and not x.startswith('.') and x in taxon_pool)]

        paralogf = os.path.join(sampledir,'all_genes_with_paralog_warnings.txt')
        paralogs = []
        with open(paralogf, 'r') as inf:
            for line in inf:
                paralogs.append(line.strip())

        print("Retrieving {} - {} (paralogs) genes from {} samples".format(len(target_genes), len(paralogs), len(sample_names)))

        gene_count = 0
        sample_count = []
        for gene in target_genes:
            if gene in paralogs:
                continue
            else:
                gene_seqs = []
                for rec in gene_seqs:
                    rec.id = rec.id.split("-")[0]
                    rec.description = ''

                for sample in sample_names:
                    if seq_dir == 'intron':
                        sample_path = os.path.join(sampledir, sample, gene, sample, 'sequences', seq_dir,"{}_{}.fasta".format(gene, filename))
                    else:
                        sample_path = os.path.join(sampledir, sample, gene, sample, 'sequences', seq_dir, gene+'.'+seq_dir)

                    if os.path.isfile(sample_path):
                        for record in SeqIO.parse(sample_path, 'fasta'):
                            record.id = convert_names[record.id[0:10]]
                            record.description = ''
                            gene_seqs.append(record)
                            sample_count.append(record.id)
                        #gene_seqs.append(SeqIO.read(sample_path, 'fasta'))

                print("Found {} sequences for {}.".format(len(gene_seqs), gene))

                if seq_dir == 'intron':
                    outfilename = "{}_{}.fasta".format(gene, filename)
                else:
                    outfilename = gene + '.' + seq_dir

                if gene_seqs:
                    SeqIO.write(gene_seqs, open(outfilename,'w'), 'fasta')
                    gene_count += 1

        if r == 'dna':
            os.system('mv *.FNA {}-hybpiper_sequences/dna'.format(STEP))
        elif r == 'intron':
            os.system('mv *.fasta {}-hybpiper_sequences/intron'.format(STEP))
        elif r == 'supercontig':
            os.system('mv *.fasta {}-hybpiper_sequences/supercontig'.format(STEP))
        elif r == 'aa':
            os.system('mv *.FAA {}-hybpiper_sequences/aa'.format(STEP))

        print('\n{} sequences recovered: {}'.format(r, gene_count))
        print('{} samples recovered: {}\n'.format(len(set(sample_count)), ','.join(list(set(sample_count)))))





def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-hy', '--hybpiper', help='path to hybpiper output directory', required=True)
    parser.add_argument('-ta', '--taxon_filtering', help='path to list of taxa to include or exclude, if None then assume to include all taxa')
    parser.add_argument('-tm', '--taxon_filtering_mode', help='\"include\" or \"exclude\" taxon filtering mode')
    #parser.add_argument('-ll', '--locus_length', action='store_true', help='include to filter by locus length')
    #parser.add_argument('-ml', '--missing_loci', action='store_true', help='include to produce trees/alignments combining missing loci (supercontig + exon)')
    parser.add_argument('-rt', '--rogue_taxa', action='store_true', help='include to filter out rogue_taxa')
    parser.add_argument('-ar', '--all_rogue_taxa', action='store_true', help='include to include ALL rogue taxa')

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()
    print(args)
    STEP = 1

    ##### RETRIEVE SEQEUENCES
    ####################################################################################################################
    taxonnames = '/projects/clement-lab/5-HybSeq/2-hybpiper_Dipsacales/Dipsacales_Baits_Project/taxonnames.csv'
    rogue_taxa_f = '/projects/clement-lab/5-HybSeq/2-hybpiper_Dipsacales/DipsacalesHybSeq/filtering_scripts/rogue_taxa.txt'
    all_taxa = []
    convert_names = {}
    names_convert = {}
    include_taxa = []
    exclude_taxa = []
    rogue_taxa = []

    with open(taxonnames, 'r') as inf:
        for line in inf:
            all_taxa.append(line.strip().split(',')[0])
            convert_names[line.strip().split(',')[0]] = line.strip().split(',')[1]
            names_convert[line.strip().split(',')[1]] = line.strip().split(',')[0]

    if args.rogue_taxa:
        with open(rogue_taxa_f, 'r') as inf:
            for line in inf:
                rogue_taxa.append(line.strip().split(',')[0])
                
        exclude_taxa.extend(rogue_taxa)

    if not args.taxon_filtering:
        pass
    else:
        with open(os.path.abspath(args.taxon_filtering), 'r') as inf:
            if args.taxon_filtering_mode == 'include':
                for line in inf:
                    taxon = line.strip()
                    if names_convert[taxon] not in include_taxa:
                        include_taxa.append(names_convert[taxon])

            elif args.taxon_filtering_mode == 'exclude':
                for line in inf:
                    taxon = line.strip()
                    if names_convert[taxon] not in exclude_taxa:
                        if names_convert[taxon] in rogue_taxa and not args.all_rogue_taxa:
                            exclude_taxa.append(names_convert[taxon])
                        elif names_convert[taxon] in rogue_taxa and args.all_rogue_taxa:
                            pass
                        elif names_convert[taxon] not in rogue_taxa:
                            exclude_taxa.append(names_convert[taxon])

    #if args.all_rogue_taxa:
    #    with open(rogue_taxa, 'r') as inf:
    #        for line in inf:
    #            if line.strip().split(',')[0] not in include_taxa:
    #                include_taxa.append(line.strip().split(',')[0])

    taxon_pool = []
    if not include_taxa:
        taxon_pool = [i for i in all_taxa if i not in exclude_taxa]
    else:
        taxon_pool = [i for i in all_taxa if i in include_taxa and i not in exclude_taxa]


    os.system('mkdir {}-hybpiper_sequences'.format(STEP))
    os.system('mkdir {}-hybpiper_sequences/dna {}-hybpiper_sequences/intron {}-hybpiper_sequences/supercontig {}-hybpiper_sequences/aa'.format(STEP, STEP, STEP, STEP))
    retrieve_sequences(os.path.abspath(args.hybpiper), taxon_pool, STEP, convert_names)

    # FIXME: generate stats for how many taxa/loci per group, how many bases total

    STEP += 1



    ##### GUIDANCE FILTERING
    ####################################################################################################################

#    os.system('mkdir {}-guidance'.format(STEP))
#    os.system('mkdir {}-guidance/dna'.format(STEP))
#    os.system('mkdir {}-guidance/supercontig'.format(STEP))
#    os.system('mkdir {}-guidance/intron'.format(STEP))
#
#    os.system()




if __name__ == '__main__':
    main()
