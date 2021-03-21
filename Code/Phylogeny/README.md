# DipsacalesHybSeq
Dipscales hybrid-sequencing project using targeted enrichment of 353 low-copy nuclear loci in the Angiosperms353 bait kit

### Methods
steps after assembly in HybPiper
1. retrieve aa/exon/intron/supercontig sequences using `retrieve_hybpiper_sequences.py`
   - this is a modified version of `retrieve_sequences.py` from the original HybPiper repo
2. align using pasta (`run_pasta.py`)
3. remove extremely gappy regions using phyutility clean (`clean_aln.py`, cutoff=0.05)
4. infer gene trees using raxml (`raxml.sh`)
5. trim spuriously long tips from trees/alignments (`trim_long_tips.py`, cutoff=0.45)
6. align cleaned/trimmed alignments using pasta (`run_pasta.py`)
7. remove gappy regions less conservatively using phyutility clean, producing final alignments (`clean_aln.py`, cutoff=0.30)
8. infer final gene trees using raxml (`raxml.sh`)
9. collapse nodes with bootstrap support below a threshold (`collapse_nodes.py`, cutoff=33)
10. combine alignments and trees into exon-only, supercontig-only, and exon+supercontig datasets (`missing_loci.py`)

RAxML tree inference:
1. concatenate alignments using AMAS
2. infer trees using `RAXML-HPC -f a ...` (``) (500 bootstraps)

ASTRAL tree inference:
1. concatenate gene trees using (`cat`)
2. infer trees using default ASTRAL settings (default: 100 bootstraps)