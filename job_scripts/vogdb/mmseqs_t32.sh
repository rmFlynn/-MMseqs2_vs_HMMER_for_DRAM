#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=10gb
#SBATCH --time=1-00:00:00
#SBATCH --job-name=vog_mmseqs
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rmflynn@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

eval "$(conda shell.bash hook)"

conda activate DRAM

cd /home/projects/DRAM/hmmer_mmseqs2_testing_take_2/shale_viral_sub_project/job_scripts/

python mmseqs.py --threads 32 \
	--working_dir "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/job_scripts/vogdb"
	--protein_data "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/vogdb/vog.proteins.all.20210525.fa" \
	--raw_algs "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/vogdb/vog.raw_algs.20210525.tar.gz" \
	--name "vogdb_t32"
