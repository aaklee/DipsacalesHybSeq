# DipsacalesHybSeq
Dipscales hybrid-sequencing project using targeted enrichment of 353 low-copy nuclear loci in the [Angiosperms353 bait kit](https://arborbiosci.com/genomics/targeted-sequencing/mybaits/mybaits-expert/mybaits-expert-angiosperms-353/) designed by [Johnson _et al._ (2019)](https://academic.oup.com/sysbio/article/68/4/594/5237557). All work in this repository is associated with 

> Lee, Aaron K., I. S. Gilman, M. Srivastav, A. Lerner, M. J. Donoghue, W. L. Clement. (*Accepted*) **Reconstructing Dipsacales phylogeny using Angiosperms353: Issues and insights**. American Journal of Botany.

Please cite us appropriately if you find our results or methodology useful!

## Using this repo
In this repository you'll find all scripts and notebooks for reproducing the results presented in Lee _et al._ (####), as well the most of the outputs (including alignments, gene and species trees, and figures). Notebooks are statistical analyses and self contained, and if they give you a "Sorry, something went wrong. Reload?" error, try copying and pasting the information in https://nbviewer.jupyter.org. Most other analyses (e.g., alignment, tree inference, and divergence time estimation) were completed on clusters at Yale University and The College of New Jersey, and can be walked through in detail in the [Wiki](https://github.com/aaklee/DipsacalesHybSeq/wiki).

## Methods
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

## References
- Matthew G Johnson, Lisa Pokorny, Steven Dodsworth, Laura R Botigué, Robyn S Cowan, Alison Devault, Wolf L Eiserhardt, Niroshini Epitawalage, Félix Forest, Jan T Kim, James H Leebens-Mack, Ilia J Leitch, Olivier Maurin, Douglas E Soltis, Pamela S Soltis, Gane Ka-shu Wong, William J Baker, Norman J Wickett, A Universal Probe Set for Targeted Sequencing of 353 Nuclear Genes from Any Flowering Plant Designed Using k-Medoids Clustering, Systematic Biology, Volume 68, Issue 4, July 2019, Pages 594–606, https://doi.org/10.1093/sysbio/syy086
