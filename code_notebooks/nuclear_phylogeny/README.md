# DipsacalesHybSeq
Dipscales hybrid-sequencing project using targeted enrichment of 353 low-copy nuclear loci in the Angiosperms353 bait kit

### Methods
steps after assembly in HybPiper
1. retrieve aa/exon/intron/supercontig sequences using `retrieve_hybpiper_sequences.py`
   - this is a modified version of `retrieve_sequences.py` from the original HybPiper repo
2. align using pasta (`run_pasta.py`)
3. remove extremely gappy regions using phyutility clean (`clean_aln.py`, cutoff=0.05)
4. infer gene trees using raxml (`raxml.sh`)
5. trim long tips from tree/alignments (`setup_treeshrink.py`, b=0.1)
6. align cleaned/trimmed alignments using pasta (`run_pasta_after_treeshrink.py`)
7. remove gappy regions less conservatively using phyutility clean, producing final alignments (`clean_aln.py`, cutoff=0.30)
8. infer final gene trees using raxml (`raxml.sh`)
9. collapse nodes with bootstrap support below a threshold (`collapse_nodes.py`, cutoff=33)
10. combine alignments and trees exon-only and exon+flanking datasets (`tree_inference_set.py`)

RAxML tree inference:
1. concatenate alignments using AMAS concat
2. add RAxML partition format to partitions file: `sed -i -e "s/^/DNA, /g" *partitions.txt`
3. infer trees using `RAXML-HPC -f a ...` (``)

ASTRAL tree inference:
1. concatenate gene trees: `cat *collapsed.gene.tre`
2. infer trees using default ASTRAL settings

Alternative uses:
- tip trimming step (5): absolute branch length values can be trimmed using `trim_long_tips.py` and test different cutoffs
- combining alignments (10): if you want to find the overlap between exon-only and exon+flanking datasets, use `missing_loci.py`
- parallel instances of pasta (`run_pasta_parallel.py`)
- parallel instances of raxml (`run_rax.py`)

Additional files:
- `taxonnames.csv` conversion from sample names to species names
- `taxa_excluded.txt` list of taxa to exclude from analyses
