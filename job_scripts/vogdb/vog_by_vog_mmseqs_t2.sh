#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=10gb
#SBATCH --time=3-00:00:00
#SBATCH --job-name=vogdb_by_vogdb_mmseqs_t2
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rmflynn@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

eval "$(conda shell.bash hook)"

conda activate DRAM

python /home/projects/DRAM/hmmer_mmseqs2_testing_take_3/src/mmseqs.py --threads 2 \
	--output_dir "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/vogdb_by_vogdb/mmseqs/" \
	--protein_data "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/vogdb/vog.proteins.all.20210525.fa" \
	--raw_aln "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/vogdb/vog.raw_algs.20210525.tar.gz" \
	--name "vogdb_by_vogdb_mmseqs_t2" \
	--clock_run True
