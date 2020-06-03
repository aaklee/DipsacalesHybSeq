#!/usr/bin/env python

#++++++++++++++++++++++++++++++ GUIDANCE filtering +++++++++++++++++++++++++++++
#
# This script 1) visualizes and analyzez patterns of alignment quality and 2)
# writes filtered fasta files based on user-defined quality thresholds.
#
# Files of interest in every GUIDANCE run:
#   MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names :Gives sequence score for
#        each sample
#   MSA.MAFFT.Guidance2_res_pair_col.scr.csv: Gives each column score
#   MSA.MAFFT.Guidance2_msa.scr: Some metadata
#
#                                                               11 December 2018
#                                                                     Ian Gilman

import re, os, sys, errno
try: from optparse import OptionParser
except ImportError:
    print("\n\tError: OptionParser (optparse) is not installed/loaded")
    sys.exit()
from glob import glob
from subprocess import Popen, PIPE
import numpy as np
import pandas as pd
from collections import Counter

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_guidance_dir(directory):
    '''Parses the relevant information in a GUIDANCE output directory.

    Parameters
    ----------
    directory <str> : path to GUIDANCE output directory'''

    seq_score_path = os.path.join(directory, "MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names")
    col_score_path = os.path.join(directory, "MSA.MAFFT.Guidance2_res_pair_col.scr.csv")
    metadata_path = os.path.join(directory, "MSA.MAFFT.Guidance2_msa.scr")
    fasta_path = os.path.join(directory, "MSA.MAFFT.aln.With_Names")

    with open(metadata_path, "r") as m:
        lines = m.readlines()
        for line in lines:
            if "#COLUMNS" in line: no_cols = int(line.split()[-1])

    # Read sequence scores and drop file ending
    seq_scores = pd.read_csv(seq_score_path, sep="\t")
    seq_scores = seq_scores.loc[seq_scores["SEQUENCE_NAME"]!="#END"]

    # Read column scores
    col_scores = pd.read_csv(col_score_path)
    col_scores.rename({"#COL_NUMBER":"column", "#RES_PAIR_COLUMN_SCORE":"score"}, axis=1, inplace=True)
    # Create analogous DF for sites with single taxon (these recieve no score from GUIDANCE)
    single_cols = pd.DataFrame(columns=col_scores.columns)
    single_cols["column"] = list(set(range(1, no_cols+1))-set(col_scores.column))
    single_cols["score"] = np.nan
    # Append single taxon sites and reset indices
    col_scores = col_scores.append(single_cols)
    col_scores.sort_values("column", inplace=True)
    col_scores.reset_index(inplace=True, drop=True)

    return seq_scores, col_scores, fasta_path



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def deinterleaf(fasta):
    '''Converts an interleaved fasta file to a non-interleaved fasta

    Parameters
    ----------
    fasta <str> : path to fasta file'''

    with open(fasta, 'r') as f:
        lines = f.readlines()

    output_data = []
    for i, line in enumerate(lines):
        if line.startswith('>'):
            sample = line
            j = i+1
            sequence = str()
            while j<len(lines) and not lines[j].startswith('>'):
                sequence+=lines[j]
                j+=1
            sequence = sequence.replace('\n', '')
            sequence = sequence.replace('\r', '')
            output_data.append(sample)
            output_data.append(sequence+'\n')

    return output_data

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def count_fasta_samples(fasta):
    '''Counts the number of samples in a fasta file

    Parameters
    ----------
    fasta <str> :  path to a fasta file'''

    with open(fasta, 'r') as f:
        lines = f.readlines()
        return len([l for l in lines if l[0] == '>'])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def remove_sequences(fasta, to_remove, as_list=False):
    '''Removes sequences from a fasta file

    Parameters
    ----------
    fasta <str> : path to fasta file
    to_remove <list-like> : a list of strings indicating sequence names to remove
    as_list <bool> : indicate if passing a list with the fasta data instead of a path to it'''

    if not as_list:
        # Deinterleaf
        flat_fasta = deinterleaf(fasta=fasta)
    else:
        flat_fasta = fasta

    # Reshape, pairing names with sequences
    rows = len(flat_fasta) / 2
    fasta_array = np.reshape(flat_fasta, (int(rows), 2))

    # Mark samples to be removed from array
    to_remove = list(to_remove)
    retained_samples = [not sample[0].replace(">", "").replace("\n", "") in to_remove for sample in fasta_array]
    # Drop samples from array
    retained_fasta_array = fasta_array[retained_samples]
    # Flatten array to list
    retained_fasta = list(np.reshape(retained_fasta_array, (np.shape(retained_fasta_array)[0]*np.shape(retained_fasta_array)[1],)))

    seqs_removed = len(fasta_array)-len(retained_fasta_array)

    if len(retained_fasta) / 2 < 4:
        print("WARNING: RESULTING FILE CONTAINS LESS THAN 4 SEQUENCES")
    if seqs_removed == 0:
        print("NO SEQUENCES REMOVED")

    return retained_fasta, seqs_removed



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def threshold_thin_sequence(sequence, scores, threshold):
    '''Removes sites from a sequence based on their guidance column
    scores and a threshold

    Parameters
    ----------
    sequence <str> : sequence string
    scores <list-like> : list of scores
    threshold <float> : column score cutoff'''
    sequence = sequence.replace("\n", "")

    return ''.join(ch for (ch, score) in zip(sequence, scores) if (score >= threshold or score==np.nan))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def remove_columns(fasta, col_scores, threshold):
    '''Removes sites from all sequences in a fasta file using
    threshold_thin_sequence

    Parameters
    ----------
    fasta <str> : path to fasta file
    col_scores <list-like> : list of column scores
    threshold <float> : column score cutoff'''

    thinned = []
    for line in deinterleaf(fasta):
        if line.startswith(">"):
            thinned.append(line)
        else:
            thinned.append(threshold_thin_sequence(line, col_scores["score"], threshold=threshold)+"\n")
    cols_removed = len(line)-len(thinned[-1])
    return thinned, cols_removed




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def main():

    print('\n~~~~~~ GUIDANCE filtering ~~~~~~\n')

    parser = OptionParser(prog="guidance_filtering", usage="%prog [options]", version="%prog 1.0")
    parser.add_option("-d", "--guidance-dir",
                      action="store",
                      dest="guidance_dir",
                      help="path GUIDANCE directory")
    parser.add_option("-c", "--column-threshold",
                      action="store",
                      type="float",
                      dest="column_threshold",
                      help="column quality threshold")
    parser.add_option("-s", "--sequence-threshold",
                      action="store",
                      type="float",
                      dest="sequence_threshold",
                      help="sequence quality threshold")
    parser.add_option("-C", "--remove-cols",
                      action="store_true",
                      dest="remove_cols",
                      help="remove columns below column threshold in output")
    parser.add_option("-S", "--remove-seqs",
                      action="store_true",
                      dest="remove_seqs",
                      help="remove sequences below sequence threshold in output")

    (options, args) = parser.parse_args()

    if not all((options.guidance_dir, options.column_threshold, options.sequence_threshold)) and (options.remove_cols or options.remove_seqs):
        print("\n\tMust specify guidance directory with -d, column threshold with -c, and sequence threshold with -s.\n")
        print("\tMust specify whether to remove columns (-C), sequences (-S), or both (-C -S).\n")
        sys.exit()


    # Parse inputs
    guidance_dir = options.guidance_dir
    column_threshold = options.column_threshold
    sequence_threshold = options.sequence_threshold

    if options.remove_cols:
        remove_cols = True
    else:
        remove_cols = False

    if options.remove_seqs:
        remove_seqs = True
    else:
        remove_seqs = False

    # Parse guidance data
    seq_scores, col_scores, fasta_path = parse_guidance_dir(guidance_dir)
    to_remove = seq_scores.loc[seq_scores["SEQUENCE_SCORE"] < sequence_threshold]["SEQUENCE_NAME"]

    # Get locus name from guidance directory
    locus = re.split("/", os.path.abspath(guidance_dir))[-1]

    if remove_seqs and remove_cols:
        # Remove columns
        thinned_cols, cols_removed = remove_columns(fasta = fasta_path, col_scores=col_scores, threshold=column_threshold)
        # Remove sequences
        passing_seqs, seqs_removed = remove_sequences(fasta=thinned_cols, to_remove=to_remove, as_list=True)
        # Write output fasta
        output_fasta_name = "{}_seq{:1.0f}_col{:1.0f}.fasta".format(locus, sequence_threshold*100,column_threshold*100)
        output_fasta_path = os.path.join(guidance_dir, output_fasta_name)
        with open(output_fasta_path, "w+") as o:
            o.writelines(passing_seqs)
        print(("Removed {} columns and {} sequences and wrote output fasta to {}".format(cols_removed, seqs_removed, output_fasta_path)))

    elif remove_seqs:
        # Remove sequences
        passing_seqs, seqs_removed = remove_sequences(fasta=fasta_path, to_remove=to_remove)
        # Write output fasta
        output_fasta_name = "{}_seq{:1.0f}.fasta".format(locus, sequence_threshold*100)
        with open(output_fasta_name, "w+") as o:
            o.writelines(passing_seqs)
        print(("Removed {} sequences and wrote output fasta to {}".format(seqs_removed, output_fasta_name)))

    elif remove_cols:
        # Remove columns
        thinned_cols, cols_removed = remove_columns(fasta = fasta_path, col_scores=col_scores, threshold=column_threshold)
        # Write output fasta
        output_fasta_name = "{}_col{:1.0f}.fasta".format(locus, column_threshold*100)
        output_fasta_path = os.path.join(guidance_dir, output_fasta_name)
        with open(output_fasta_path, "w+") as o:
            o.writelines(thinned_cols)
        print(("Removed {} columns and wrote output fasta to {}".format(cols_removed, output_fasta_path)))

    else:
        print("No sequences or columns removed. To remove, pass flags --remove_seqs and/or --remove_cols")



#+------------------------------------------------------------------------
if __name__ == "__main__":
    main()
