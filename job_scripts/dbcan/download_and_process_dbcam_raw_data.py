"download and process the bdcam db using dram tools"
import os

os.chdir("/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/dbcan/")

# Manually get hmmer
os.system("wget http://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V9.txt")
os.system("hmmpress -f dbCAN-HMMdb-V9.txt")

# get additional files we may need
os.system("wget http://bcb.unl.edu/dbCAN2/download/dbCAN-fam-aln-V9.tar.gz")
os.system("wget http://bcb.unl.edu/dbCAN2/download/CAZyDB.07312020.fa")

