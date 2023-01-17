import subprocess
import argparse
import sys

parser = argparse.ArgumentParser(description="Running the postprocess")
parser.add_argument('data_dir', type=str, help='Input dir which postprocess folder located')
args = parser.parse_args()
path_to_dir = args.data_dir+"/postprocess/"
sys.path.append(path_to_dir)
from post_process_input import *

print('')
print("<----- Make Cipher Files ----->")
subprocess.call(['python',path_to_dir+'make_ciphers.py'])
print('')
print('')
print("<----- Make Result Files ----->")
subprocess.call(['python',path_to_dir+'make_results.py'])
print('')
print('')
print("<----- Make Report Files ----->")
subprocess.call(['python',path_to_dir+'make_reports.py'])

