#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=10gb
#SBATCH --time=1-00:00:00
#SBATCH --job-name=cazydb_by_dbcan_mmseqs
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rmflynn@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

eval "$(conda shell.bash hook)"

conda activate DRAM

python /home/projects/DRAM/hmmer_mmseqs2_testing_take_3/src/mmseqs.py --threads 32 \
	--output_dir "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/cazydb_by_dbcan/mmseqs/" \
	--protein_data "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/dbcan/CAZyDB.07312020.fa" \
	--raw_aln "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/dbcan/dbCAN-fam-aln-V9.tar.gz" \
	--name "cazydb_by_dbcan_mmseqs" \
	--clock_run False
