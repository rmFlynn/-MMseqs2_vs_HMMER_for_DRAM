#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=8gb
#SBATCH --time=240:00:00
#SBATCH --job-name=genome15_by_vogdb_hmmer
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rmflynn@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

eval "$(conda shell.bash hook)"

conda activate DRAM

python /home/projects/DRAM/hmmer_mmseqs2_testing_take_3/src/hmmer.py \
 	--output_dir "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/genome15_by_vogdb/hmmer/" \
 	--protein_data "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/genome15/genes.faa" \
	--hmm "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/vogdb/vog.hmm/vog.hmm.20210525.txt" \
 	--name 'genome15_by_vogdb_hmmer' \
 	--clock_run Flase
