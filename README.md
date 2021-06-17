# HMMER3 VS MMseqs2 Testing and Evaluation Project

This folder is yet another reorganization of the HMMER and MMseqs testing project. The goal is evaluating the value of replacing HMMER in DRAM with MMseqs to improve the speed and flexibly of DRAM.

The organization of this folder is guided by the fallowing principles.

 * One truth: There can be only version of any data file, script, or function that has the same name or clames the same goal.
 * Don't original: All files in the data folder are original and never modified, there source is documented bellow
 * Git: All code/scrips and changes must be in git FULL STOP

These rules are simplifications/bastardizations of sage advice / things that have been yelled at me.

There is further reading:
  * https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005510

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

For more information on organization see the log.txt


## log.txt

This will list all test with the version of code and scripts used to run them.
