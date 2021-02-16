#!/usr/bin/bash

#SBATCH --workdir=./
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH --job-name=raxml
#SBATCH --output=raxml_%j.log
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --partition=long
#SBATCH --time=2-00:00:00


# get a list of alignment files
FILES=$(ls /projects/clement-lab/5-HybSeq/2-hybpiper_Dipsacales/H-clustering_first2_selective/7-phyutility/phyutility30/0.45/supercontig05/*.fasta)

let COUNT=1
for i in $FILES;
do
    # parse alignment name
    # first cut: separate by "/" and take last entry
    # second cut: take first four characters of fasta file name (locus id)
    j=$(ls $i | cut -d/ -f11 | cut -c1-4)

    # write to commands.txt file
    echo "module add raxml; raxmlHPC -fa -x 65135 -# 100 -p 65135 -m GTRGAMMA -s $i -n $j.rax$COUNT"
    let COUNT++

done > commands.txt

# run in parallel using GNU parallel
parallel < commands.txt
