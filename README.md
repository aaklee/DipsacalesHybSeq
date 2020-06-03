# DipsacalesHybSeq
Dipscales hybrid-sequencing project

### Methods
steps after assembly in HybPiper
1. retrieve aa/exon/intron/supercontig sequences using `retrieve_hybpiper_sequences.py`
   - this is a modified version of `retrieve_sequences.py` from the original HybPiper repo
2. align using MAFFT, score using guidance (guidance will align and score)
3. filtering sequences and columns in alignments using `guidance_filtering-AL.py3`
4. concatenate individual locus alignments into a supermatrix using AMAS
5. reconstruct individual locus phylogenies using RAxML
6. collapse nodes with bootstrap < 33 using `collapse_nodes.py`
7. reconstruct species phylogenies under concatenated (RAxML) or coalescent (ASTRAL-III) approach
