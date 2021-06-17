#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=8gb
#SBATCH --time=240:00:00
#SBATCH --job-name=shale_hmmer
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rmflynn@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low


eval "$(conda shell.bash hook)"

conda activate DRAM

cd /home/projects/DRAM/hmmer_mmseqs2_testing_take_2/shale_viral_sub_project/job_scripts/

	python /home/projects/DRAM/hmmer_mmseqs2_testing_take_3/src/hmmer.py \
		--hmm "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/vogdb/vog.hmm/vog.hmm.20210525.txt" \
         	--genes_faa "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/vogdb/vog.proteins.all.20210525.fa" \
	 	--output_dir "/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/results/vogdb/hmmer" \
	 	--name 'vog_hmmer' \
	 	--clock_run True
