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

  a. Pipelines/Sub-project

2. Results: Time stamped results

  a. Pipelines/Sub-project

3. Analysis: Analysis of results

  a. Pipelines/Sub-project
  b. Shared: code that is used across sub projects

4. Job_scripts: Writon-ization of a scripts file for slurm jobs and scripts for jobs

  a. Pipelines/Sub-project

5. src: Code referenced by job scripts, mostly python.

  a. Sub-project
  b. Shared: code that is used across sub projects

## Data sets (location:Origin)

   * dbcam/
       * CAZyDB.07312020.fa : http://bcb.unl.edu/dbCAN2/download/dbCAN-fam-aln-V9.tar.glz
       * fadbCAN-fam-aln-V9.tar.gz : http://bcb.unl.edu/dbCAN2/download/CAZyDB.07312020.fa
       * dbCAN-HMMdb-V9.txt : http://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V9.txt
   * genome15/
        * genes.faa : /home/projects-wrighton/DRAM_performance/15_soil_genomes/DRAM_test/genes.faa
   * shale/
        * genes.faa : /home/projects/DRAM/Shale_10kb_pub_050420/DRAMv_nov/genes.faa
   * vogdb/ : http://fileshare.csb.univie.ac.at/vog/vog203/
        * vog.hmm/ : vog.hmm.20210525.tar.gz : http://fileshare.csb.univie.ac.at/vog/vog203/vog.hmm.tar.gz
        * vog.members.20210525.tsv.gz : http://fileshare.csb.univie.ac.at/vog/vog203/vog.members.tsv.gz
        * vog.proteins.all.20210525.fa : http://fileshare.csb.univie.ac.at/vog/vog203/vog.proteins.all.fa.gz
        * vog.raw_algs.20210525.tar.gz : http://fileshare.csb.univie.ac.at/vog/vog203/vog.raw_algs.tar.gz

Results

The results folder is organized first by the pipeline(aka subproject) then by MMseqs or HMMER.
The results files are also named and dated to avoid confusion.
Note that the pipelines are cazydb_by_dbcan, genome15_by_dbcan, genome15_by_vogdb, shale_by_vogdb, and vogdb_by_vogdb.

## Analysis
The analysis is the most complex part of the pipeline, it has the most code and so the most opportunity for things to go wrong.
In the analysis folder there are folders for code(the code used for analysis), data(data tables made by analysis), and output(tables and figures that are used in the presentation. Each of these folders contains folders for each pipeline.

### Code
Everything here hinges on the data in the results folder being read and processed correctly and so an effort was made to check each pipeline for errors in the ID fields.”RECHECK IN PROGRESS”
To look at these tests open the test_id_validity.py folder in each pipeline folder.
Only one MMSeqs model was tested, specifically the e-value = 1e-15, sensitivity=7.5 compared to a HMMER model with the same e-value. Both models were filtered by percent coverage.

In the pipelines that were based off of VOGdb (genome15_by_vogdb, shale_by_vogdb, and vogdb_by_vogdb)  the following was tested:
Protein ids appear only once.
No protein ids in the output of mmseqs or hmmer were not from the input faa file. The vog truth was also checked this way.
No Dogs that were not in the raw alignments were in any outputs.
What input protein ids or vog ids appeared in no matches.
The only broad note from this analis was that mmseqs left out slightly more proteins and vogs compared to HMMER.
In the dbCAN based pipelines (cazydb_by_dbcan and genome15_by_dbcan) the following was tested:
Protein ids appear more than once. Unlike in vog.
No protein ids in the output of mmseqs or hmmer were not from the input faa file. Like in vog.
What input protein ids had no matches. Like in vog.
No bANds that were not in the raw alignments were in any outputs. Like in vog.
What ids were included and excluded by the regex filter.
Overall takeaways from this analysis are: MMseqs was missing many more proteins from the input in the output meaning they were not matched. The regex filter removed some valid proteins that were not in dbCAN. This filter was only used on cazydb_by_dbcan so proteins outside the scope of dbCAN were dropped, and CBM35inCE17 was also dropped because it was not clear how it would be matched.

#### Output
This folder contains carefully organized plots and tables used in the presentation. It should be self explanatory but more information may be on the way.
#### Data
Along with plots and tables for the presentation there are also tables of summary statistics. These are used in some of the tables and should be run first in some pipelines. There are two types of output tables, gold standard, and hmmer vs mmseqs comparisons. The hmmer vs MMseqs comparisons are in all pipelines and simply count how many times mmseqs agreed or disagreed with the model with each set of settings. The Gold Standard stats are only in  cazydb_by_dbcan and vogdb_by_vogdb, and are unique to each case in what they included. Both Gold Standard stats tables compare known true andotations to the predicted annotations of the models.

Both folders contain all data in both csv and pkl format.

### Pipeline Results Processing Summary

All the pipelines are different, though they share similarities based on whether they are andotated by VOGdb or dbCAN. The most critical part of the results analisis is the post processing that changes how the results are read and interpreted. Below is a summary of how the reading is different for each pipeline.
  * Genome15_by_vogdb:
     * For MMseqs:
          Remove '.msa' from all annotations
          Take only one value for each protein id, the value with the lowest e-value
     * For HMMER:
          Take only one annotations for each protein id, the value with the lowest e-value

  * Shale_by_vogdb:

     * For MMseqs:
          Remove '.msa' from all annotations
          Take only one value for each protein id, the value with the lowest e-value
     * For HMMER:
          Take only one annotations for each protein id, the value with the lowest e-value

  * Vogdb_by_vogdb:

     * For MMseqs:
          Remove '.msa' from all annotations
          Take only one value for each protein id, the value with the lowest e-value
     * For HMMER:
          Take only one annotations for each protein id, the value with the lowest e-value


   * Cazydb_by_dbcan:

     * For MMseqs:
          Remove '.aln' from all annotations
          Take multiple annotations for each protein id, but still take only the pair(protein and annotation) with the lowest e-value
          Split any annotations in the form DDD_NN to take only the portion before the under score(DDD).
          compare all resulting IDs to a filtering regex and take only matches
     * For HMMER:
          Take multiple annotations for each protein id, but still take only the pair(protein and annotation) with the lowest e-value
          Remove '.hmm' from all annotations
          Split any annotations in the form DDD_NN to take only the portion before the under score(DDD).
          compare all resulting IDs to a filtering regex and take only matches

  * Genome15_by_dbcan:

     * For MMseqs:
          Remove '.aln' from all annotations
          Take multiple annotations for each protein id, but still take only the pair(protein and annotation) with the lowest e-value
          Split any annotations in the form DDD_NN to take only the portion before the under score(DDD).
          DON'T compare resulting IDs to a filtering regex, there is no gold truth in this case only HMMER VS MMseqs comparison
     * For HMMER:
          Take multiple annotations for each protein id, but still take only the pair(protein and annotation) with the lowest e-value
          Remove '.hmm' from all annotations
          Split any annotations in the form DDD_NN to take only the portion before the under score(DDD).
          DON'T compare resulting IDs to a filtering regex, there is no gold truth in this case only HMMER VS MMseqs comparison


## TO DO

Finish all major tests.

Validate data choussess for reading and filtering data.

Add additional info to readme

Finish the time plots

Clear all plots and analysis data and re-run

