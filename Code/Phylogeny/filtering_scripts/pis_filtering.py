#!/usr/bin/python

# filtering of LOCI based on the proportion of parsimony-informative sites (PIS) in an alignment
# can be run on pre- or post-alignment filtering in GUIDANCE2, TRIMAL, ETC...

# IF YOU'RE RUNNING INTO ERRORS, maybe check this out:
# uses AMAS to calculate PIS counts/proportions, so make sure that you have AMAS.py in the path

# output will be a concatenated set of trees corresponding to loci with a certain proportion of PIS
# thresholds determined by the user

import sys, os
import pandas
from matplotlib import pyplot


aln_dir = os.path.abspath(sys.argv[1])
lower = float(sys.argv[2]) # lower boundary, as a decimal
upper = float(sys.argv[3]) # upper boundary, as a decimal


aln_f = []
for f in os.listdir(aln_dir):
    if f.endswith('fasta') and 'bootstrap' not in f and 'renamed' not in f:
        aln_f.append(os.path.join(aln_dir, f))


os.system('AMAS.py summary -f fasta -d dna -i {} -o AMAS_summary.out'.format(' '.join(aln_f)))


amas = pandas.read_csv('AMAS_summary.out', sep='\t', header=0)
amas['Alignment_name'] = amas['Alignment_name'].apply(lambda x: x[0:4])
aln_length = amas['Alignment_length'].sum()

print('\ntotal alignment length of {}:\t{}\n'.format(f.split('/')[-1], aln_length))

# generate histogram describing distribution of percent PIS
pis_perc = amas[['Alignment_name', 'Proportion_parsimony_informative']]

pis_hist, ax = pyplot.subplots()
ax.hist(pis_perc['Proportion_parsimony_informative'], bins=50)
ax.set_xlabel('Proportion_parsimony_informative')
ax.set_ylabel('frequency')
pis_hist.savefig('pis_histogram.pdf')


# sort pis_percentage values

pis_perc = pis_perc.sort_values(by=['Proportion_parsimony_informative'])
min_perc = pis_perc.iloc[0][1]
max_perc = pis_perc.iloc[len(pis_perc.index)-1][1]
avg_perc = pis_perc['Proportion_parsimony_informative'].sum() / len(pis_perc.index)
std_perc = pis_perc['Proportion_parsimony_informative'].std()
#med_perc = pis_perc[10][1]
print(min_perc, max_perc)
print(avg_perc)

print()


# get paralog information to avoid adding paralogs to final tree
paralogf = os.path.abspath('/projects/clement-lab/5-HybSeq/2-hybpiper_Dipsacales/rerun/1-hybpiper/trim/all_genes_with_paralog_warnings.txt')
paralogs = []
with open(paralogf, 'r') as inf:
    for line in inf:
        paralogs.append(int(line.strip()))


# actual filtering
if upper == 0.0 and lower == 0.0:
    print('no filter')
    amas_filtered = amas
elif lower > 0.0 and upper == 0.0:
    print('lower filter, no upper filter')
    amas_filtered = amas[amas['Proportion_parsimony_informative'] > lower]
elif lower == 0.0 and upper > 0.0:
    print('no lower filter, upper filter')
    amas_filtered = amas[amas['Proportion_parsimony_informative'] < upper]
else:
    print('lower and upper filters')
    amas_filtered = amas[amas['Proportion_parsimony_informative'] > lower]
    print(len(amas_filtered['Proportion_parsimony_informative']))
    amas_filtered = amas_filtered[amas_filtered['Proportion_parsimony_informative'] < upper]
    print(len(amas_filtered['Proportion_parsimony_informative']))

print('{} of {} loci remaining.\n'.format(len(list(amas_filtered['Alignment_name'])), len(list(amas['Alignment_name']))))
aln_length = amas_filtered['Alignment_length'].sum()
print('\ntotal alignment length of {}:\t{}\n'.format(f.split('/')[-1], aln_length))

# get trees to concatenate for ASTRAL analysis
to_cat = []
loci = list(amas_filtered['Alignment_name'].astype(int))
print(loci)
for f in os.listdir(aln_dir):
    if ('bestTree' in f) and ('bootstrap' not in f) and ('renamed' not in f) and (int(f.split('.')[1]) not in paralogs) and (int(f.split('.')[1]) in loci):
        to_cat.append(os.path.join(aln_dir, f))

os.system('cat {} > pis_filtering_loci.gene.tre'.format(' '.join(to_cat)))




# filtering based on stdev
print('determine filter based on stdev - get only loci corresponding to std_up and std_down')
std_up = avg_perc + std_perc
std_down = avg_perc - std_perc

print(std_down, std_up)

amas_filtered = amas[amas['Proportion_parsimony_informative'] > std_down]
amas_filtered = amas_filtered[amas_filtered['Proportion_parsimony_informative'] < std_up]

to_cat = []
loci = list(amas_filtered['Alignment_name'].astype(int))
print(loci)
for f in os.listdir(aln_dir):
    if ('bestTree' in f) and ('bootstrap' not in f) and ('renamed' not in f) and (int(f.split('.')[1]) not in paralogs) and(int(f.split('.')[1]) in loci):
        to_cat.append(os.path.join(aln_dir, f))

os.system('cat {} > stdev_pis_filtering_loci.gene.tre'.format(' '.join(to_cat)))

aln_length = amas_filtered['Alignment_length'].sum()
print('\ntotal alignment length of {}:\t{}\n'.format(f.split('/')[-1], aln_length))
