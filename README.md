# HMMER3 VS MMseqs2 Testing and Evaluation Project

This folder is yet another reorganization of the HMMER and MMseqs testing project. The goal is evaluating the value of replacing HMMER in DRAM with MMseqs to improve the speed and flexibly of DRAM.

The organization of this folder is guided by the fallowing principles.

 * One truth: There can be only version of any data file, script, or function that has the same name or clames the same goal.
 * Don't original: All files in the data folder are original and never modified, there source is documented bellow
 * Git: All code/scrips and changes must be in git FULL STOP

These rules are simplifications/bastardizations of sage advice / things that have been yelled at me.

There is further reading:
  * https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005510

## log.txt

This will list all test with the version of code and scripts used to run them.

## Organization

1. Data: Unmodified original data
  a. Sub-project
2. Results: Time stamped results
  a. Sub-project
3. Analysis: Analysis of results
  a. Sub-project
  b. Shared: code that is used across sub projects
4. Job_scripts: Writon-ization of a scripts file for slurm jobs and scripts for jobs
  a. Sub-project
5. src: Code referenced by job scripts, mostly python.
  a. Sub-project
  b. Shared: code that is used across sub projects

## Data sets (location:Origin)

data/
    dbcam/
        CAZyDB.07312020.fa : http://bcb.unl.edu/dbCAN2/download/dbCAN-fam-aln-V9.tar.gz
        fadbCAN-fam-aln-V9.tar.gz : http://bcb.unl.edu/dbCAN2/download/CAZyDB.07312020.fa
        dbcan_hmm* : download_and_process_dbcan(dbcan_hmm=None, output_dir='./', dbcan_release='9', verbose=True)
    genome15/
	genes.faa : /home/projects-wrighton/DRAM_performance/15_soil_genomes/DRAM_test/genes.faa
    shale/
        genes.faa : /home/projects/DRAM/Shale_10kb_pub_050420/DRAMv_nov/genes.faa

    vogdb/ : http://fileshare.csb.univie.ac.at/vog/vog203/
         vog.hmm/ : vog.hmm.20210525.tar.gz : http://fileshare.csb.univie.ac.at/vog/vog203/vog.hmm.tar.gz
         vog.members.20210525.tsv.gz : http://fileshare.csb.univie.ac.at/vog/vog203/vog.members.tsv.gz
         vog.proteins.all.20210525.fa : http://fileshare.csb.univie.ac.at/vog/vog203/vog.proteins.all.fa.gz
         vog.raw_algs.20210525.tar.gz : http://fileshare.csb.univie.ac.at/vog/vog203/vog.raw_algs.tar.gz


