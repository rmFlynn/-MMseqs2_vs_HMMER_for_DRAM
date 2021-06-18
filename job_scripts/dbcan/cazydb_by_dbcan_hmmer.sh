#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=8gb
#SBATCH --time=240:00:00
#SBATCH --job-name=cazydb_by_dbcan_hmmer
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rmflynn@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

eval "$(conda shell.bash hook)"

conda activate DRAM

python /home/projects/DRAM/hmmer_mmseqs2_testing_take_3/src/hmmer.py \
 	--output_dir "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/cazydb_by_dbcan/hmmer/" \
 	--protein_data "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/dbcan/CAZyDB.07312020.fa" \
	--hmm "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/dbcan/dbCAN-HMMdb-V9.txt" \
 	--name 'cazydb_by_dbcan_hmmer' \
 	--clock_run Flase
