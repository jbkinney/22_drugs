import os, sys
import subprocess
import time
# Import source files
sys.path.append('./src')
from src import metadata
from metadata import *
from src import utils

start_time = time.time()


import argparse
parser = argparse.ArgumentParser(description='Running Pipleline')
parser.add_argument('data_dir', type=str, help='Input dir for Metadata and params.py files')
args = parser.parse_args()

# Import parameters from params.py in projectx directory
sys.path.append(args.data_dir)
from params import *

# Start Global timing
start_time = time.time()
# Load metadata file
metadata_class = metadata(args.data_dir+'/'+metadata_file)
metadata_class.file_existence()
metadata_class.load_metadata()
metadata_class.validate_metadata()

# Make directories in projects
if clean_intermediate:
    utils.clean_dir(args.data_dir+'/'+interm_dir)

if clean_tmp:
    utils.clean_dir(args.data_dir+'/'+tmp_dir)

# Split fastq files
path = os.getcwd()
## Run spilt_fastq script
print('  ')
print('-> Split files         <-')
print('-> ---------------- <-')
print('-> ---------------- <-')
print('-> ---------------- <-')
if make_split_files:
    subprocess.call(['python',path+'/src/split_fastq.py', args.data_dir])

# Parse Features
print('-> Make features files <-')
print('-> ---------------- <-')
print('-> ---------------- <-')
print('-> ---------------- <-')
if make_features_files:
    subprocess.call(['python',path+'/src/parse_features.py', args.data_dir])

# Count Features
print('-> Make count files    <-')
print('-> ---------------- <-')
print('-> ---------------- <-')
print('-> ---------------- <-')
print('-> ---------------- <-')
if make_counts_files:
    subprocess.call(['python',path+'/src/count_features.py', args.data_dir])

print('-> Merge feature files <-')
print('-> ---------------- <-')
print('-> ---------------- <-')
print('-> ---------------- <-')
# Merge feature files based on their LIDs:
if merge_feature:
    subprocess.call(['python',path+'/src/merge_features.py', args.data_dir])


print('Pipeline totally took {:2.2f} seconds to be executed'.format(time.time()-start_time))
