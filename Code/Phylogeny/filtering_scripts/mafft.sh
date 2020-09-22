#!/usr/bin/bash

#SBATCH --workdir=./
#SBATCH --mail-user=leea33@tcnj.edu
#SBaATCH --mail-type=ALL
#SBATCH --job-name=mafft
#SBATCH --output=mafft_%j.log
#SBATCH --ntasks=32
#SBATCH --partition=long
#SBATCH --time=3-00:00:00

echo "START: " `date`

module add mafft

for i in /projects/clement-lab/5-HybSeq/2-hybpiper_Dipsacales/G1-clustering_first/1-hybpiper_sequences/aa/*.FAA;
do
    j=$(echo $i | cut -d/ -f9 | cut -c1-4)
    mafft $i > $j.fasta
done
echo "END: " `date`
