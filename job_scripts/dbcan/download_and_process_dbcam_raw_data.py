"download and process the bdcam db using dram tools"
from mag_annotator.database_processing import download_and_process_dbcan
import os

os.chdir("/home/projects/DRAM/hmmer_mmseqs2_testing_take_3/data/dbcan/")
download_and_process_dbcan(dbcan_hmm=None, output_dir='./',
                           dbcan_release='9', verbose=True)

# get additional files we may need
os.system("wget http://bcb.unl.edu/dbCAN2/download/dbCAN-fam-aln-V9.tar.gz")
os.system("wget http://bcb.unl.edu/dbCAN2/download/CAZyDB.07312020.fa")

