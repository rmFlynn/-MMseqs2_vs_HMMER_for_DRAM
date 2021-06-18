#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=10gb
#SBATCH --time=1-00:00:00
#SBATCH --job-name=shale_by_vogdb_mmseqs
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rmflynn@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

eval "$(conda shell.bash hook)"

conda activate DRAM

python /home/projects/DRAM/hmmer_mmseqs2_testing_take_3/src/mmseqs.py --threads 32 \
	--output_dir "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/shale_by_vogdb/mmseqs" \
	--protein_data "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/shale/genes.faa" \
	--raw_aln "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/vogdb/vog.raw_algs.20210525.tar.gz" \
	--name "shale_by_vogdb_mmseqs" \
	--clock_run False
