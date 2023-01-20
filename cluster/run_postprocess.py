import subprocess
import argparse
import sys

parser = argparse.ArgumentParser(description="Running the postprocess")
parser.add_argument('data_dir', type=str, help='Input dir which postprocess folder located')
args = parser.parse_args()
path_to_dir = args.data_dir+"/postprocess/"
sys.path.append(path_to_dir)
from post_process_input import *

#print("\n<----- Make Cipher Files ----->\n")
#subprocess.call(['python',path_to_dir+'make_ciphers.py'])
print("\n<----- Make Result Files ----->\n")
subprocess.call(['python',path_to_dir+'make_results.py'])
#print("\n<----- Make Report Files ----->\n")
#subprocess.call(['python',path_to_dir+'make_reports.py'])
#print("\n<----- Make PSI Files ----->\n")
#subprocess.call(['python',path_to_dir+'make_psi.py'])
