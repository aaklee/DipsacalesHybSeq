#!/usr/bin/python

# 12-may-2020, aaron lee
# pipeline framework to filter hybpiper results

# input: path to directory containing hybpiper output
# output: everything you ever wanted

import sys
import os
import argparse
from Bio import SeqIO

def retrieve_sequences(sampledir, target_file, taxon_pool, convert_names, exclude_paralogs): # shamelessly copied & modified from hybpiper/retrieve_sequences.py
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

        # Use gene names parsed from a bait file.
        baitfile = target_file
        target_genes_dict = SeqIO.to_dict(SeqIO.parse(baitfile,'fasta'))
        target_genes = list(set([x.id.split('-')[-1] for x in SeqIO.parse(baitfile,'fasta')]))

        sample_names = [x for x in os.listdir(sampledir) if (os.path.isdir(os.path.join(sampledir,x)) and not x.startswith('.') and x in taxon_pool)]

        # check if exclude paralogs
        paralogs = []
        if exclude_paralogs:
            paralogf = os.path.join(sampledir,'all_genes_with_paralog_warnings.txt')
            with open(paralogf, 'r') as inf:
                for line in inf:
                    paralogs.append(line.strip())
            print("Retrieving {} - {} (paralogs) genes from {} samples".format(len(target_genes), len(paralogs), len(sample_names)))
        else:
            print("Retrieving {} genes from {} samples".format(len(target_genes), len(sample_names)))


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
                            if convert_names:
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

        if gene_seqs:
            if r == 'dna':
                os.system('mkdir 1-hybpiper_sequences/dna')
                os.system('mv *.FNA 1-hybpiper_sequences/dna')
            elif r == 'intron':
                os.system('mkdir 1-hybpiper_sequences/intron')
                os.system('mv *.fasta 1-hybpiper_sequences/intron')
            elif r == 'supercontig':
                os.system('mkdir 1-hybpiper_sequences/supercontig')
                os.system('mv *.fasta 1-hybpiper_sequences/supercontig')
            elif r == 'aa':
                os.system('mkdir 1-hybpiper_sequences/aa')
                os.system('mv *.FAA 1-hybpiper_sequences/aa')

        print('\n{} sequences recovered: {}'.format(r, gene_count))
        print('{} samples recovered: {}\n'.format(len(set(sample_count)), ','.join(list(set(sample_count)))))





def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-hy', '--hybpiper', help='path to hybpiper output directory', required=True)
    parser.add_argument('-tf', '--target_file', help='path to hybpiper target sequences', required=True)
    parser.add_argument('-ta', '--taxon_filtering', help='path to list of taxa to include or exclude, if None then assume to include all taxa')
    parser.add_argument('-tm', '--taxon_filtering_mode', help='\"include\" or \"exclude\" taxon filtering mode')
    #parser.add_argument('-ll', '--locus_length', action='store_true', help='include to filter by locus length')
    #parser.add_argument('-ml', '--missing_loci', action='store_true', help='include to produce trees/alignments combining missing loci (supercontig + exon)')
    parser.add_argument('-rt', '--rogue_taxa', action='store_true', help='exclude ALL rogue_taxa')
    parser.add_argument('-ar', '--all_rogue_taxa', action='store_true', help='include ALL rogue taxa')
    parser.add_argument('-sr', '--some_rogue_taxa', action='store_true', help='include SOME rogue taxa')
    parser.add_argument('-co', '--convert', action='store_true', help='convert sample names from \"Gilman\" code')
    parser.add_argument('-ep', '--exclude_paralogs', action='store_true', help='exclude potential paralogs flagged by HybPiper')
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()
    print(args)

    ##### RETRIEVE SEQEUENCES
    ####################################################################################################################
    taxonnames = '/projects/clement-lab/resources/software/DipsacalesHybSeq/filtering_scripts/taxonnames.csv'
    rogue_taxa_f = '/projects/clement-lab/resources/software/DipsacalesHybSeq/filtering_scripts/rogue_taxa.txt'
    all_taxa = []
    convert_names = {}
    names_convert = {}
    include_taxa = []
    exclude_taxa = []
    rogue_taxa = []

    # if convert, get conversions
    if args.convert:
        with open(taxonnames, 'r') as inf:
            for line in inf:
                all_taxa.append(line.strip().split(',')[0])
                convert_names[line.strip().split(',')[0]] = line.strip().split(',')[1]
                names_convert[line.strip().split(',')[1]] = line.strip().split(',')[0]
    else:
        with open(os.path.join(os.path.abspath(args.hybpiper), 'reads_list.txt'), 'r') as inf:
            for line in inf:
                all_taxa.append(line.strip().split(' ')[0].split('/')[-1].split('_')[0])

    # if rogue taxa, get rogue taxa to exclude
    if args.rogue_taxa:
        with open(rogue_taxa_f, 'r') as inf:
            for line in inf:
                rogue_taxa.append(line.strip().split(',')[0])
        exclude_taxa.extend(rogue_taxa)

    # do not filter taxa
    if not args.taxon_filtering:
        samples = []
        with open(os.path.join(os.path.abspath(args.hybpiper), 'reads_list.txt'), 'r') as inf:
            for line in inf:
                samples.append(line.strip().split(' ')[0].split('/')[-1].split('_')[0])

        include_taxa = [i for i in samples]


    # filter taxa
    else:
        with open(os.path.abspath(args.taxon_filtering), 'r') as inf:
            if args.taxon_filtering_mode == 'include':
                for line in inf:
                    taxon = line.strip()
                    if args.convert and names_convert[taxon] not in include_taxa:
                        include_taxa.append(names_convert[taxon])

                    elif not args.convert and taxon not in include_taxa:
                        include_taxa.append(taxon)

            elif args.taxon_filtering_mode == 'exclude':
                for line in inf:
                    taxon = line.strip()
                    if args.convert and names_convert[taxon] not in exclude_taxa:
                        if names_convert[taxon] in rogue_taxa and not args.all_rogue_taxa:
                            exclude_taxa.append(names_convert[taxon])
                        elif names_convert[taxon] in rogue_taxa and args.all_rogue_taxa:
                            pass
                        elif names_convert[taxon] not in rogue_taxa:
                            exclude_taxa.append(names_convert[taxon])

                    elif not args.convert and taxon not in exclude_taxa:
                        if taxon in rogue_taxa and not args.all_rogue_taxa:
                            exclude_taxa.append(taxon)
                        elif taxon in rogue_taxa and args.all_rogueTaxa:
                            pass
                        elif taxon not in rogue_taxa:
                            exclude_taxa.append(taxon)

    if args.all_rogue_taxa:
        with open(rogue_taxa_f, 'r') as inf:
            for line in inf:
                if line.strip().split(',')[0] not in include_taxa:
                    include_taxa.append(line.strip().split(',')[0])

    if args.some_rogue_taxa:
        with open(rogue_taxa_f, 'r') as inf:
            for line in inf:
                if line.strip().split(',')[0] not in include_taxa and line.strip().split(',')[0] not in exclude_taxa:
                    include_taxa.append(line.strip().split(',')[0])


    taxon_pool = []
    if not include_taxa and exclude_taxa:
        taxon_pool = [i for i in all_taxa if i not in exclude_taxa]
    elif include_taxa and not exclude_taxa:
        taxon_pool = [i for i in all_taxa if i in include_taxa]
    elif include_taxa and exclude_taxa:
        taxon_pool = [i for i in all_taxa if i in include_taxa or i not in exclude_taxa]

    os.system('mkdir 1-hybpiper_sequences')
    retrieve_sequences(os.path.abspath(args.hybpiper), os.path.abspath(args.target_file), taxon_pool, convert_names, args.exclude_paralogs)



if __name__ == '__main__':
    main()