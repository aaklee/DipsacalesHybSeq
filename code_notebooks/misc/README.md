# scripts used for one-off analyses

## extract genic/intergenic regions from plastid alignments
- you will need:
  - reference-assembled plastid genomes aligned with a reference genome
  - annotations extracted in a CSV format from Geneious
  - annotations should include "IR", "LSC", and "SSC"

### to get genic regions from CP alignments:
1. `extract_genic_regions.py` into your directory of choice
2. `tally_regions.py` for extracted regions
3. `parse_genic_by_locus.py` with the outputs from (2)

### to get intergenic regions from CP alignments:
1. `extract_intergenic_regions.py` into your directory of choice
2. `tally_regions.py` for extracted regions
3. `parse_intergenic_by_locus.py` with the outputs from (2)

## for dating
subset loci assembled using `get_genes_nuclear.py` or `get_genes_plastid.py`