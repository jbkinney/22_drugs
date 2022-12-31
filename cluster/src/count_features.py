import sys, os, subprocess
from metadata import *
import utils

data_dir = sys.argv[1]
sys.path.append(data_dir)
# Import parameters from params.py in projectx directory
from params import *
# Load metadata file
metadata_class = metadata(data_dir+'/'+metadata_file)
metadata_class.load_metadata()
sample_sheet= metadata_class.excel_file['samples']
LIDs_sheet = metadata_class.excel_file['LIDs']
LIDs_sheet['USE'] = LIDs_sheet['USE'].astype(str)
LIDs_sheet = LIDs_sheet[LIDs_sheet['USE']=='True']
LID_list   =  LIDs_sheet['LID'].tolist()

_l = []
for _z in LID_list:
    _l.append(sample_sheet.index[sample_sheet['LID']==_z].tolist())
# Flatten the list
idx=sum(_l,[])

sample_list = sample_sheet['sample'][idx]

path = os.getcwd()

commands = [] 
# For each sample collect the counts from the feature files
for sample_name in sample_list:
    # The following is the python command to make 
    # count files for each sample
    #subprocess.call(['python',
    #                  path+'/src/make_counts_files.py',
    #                  '-d', path+'/'+data_dir,
    #                  '-s', str(sample_name)])

    command = 'python %s -d %s -s %s'% \
              (path+'/src/make_counts_files.py',
               path+'/'+data_dir, 
               str(sample_name))
    commands.append(command)

utils.submit_and_complete_jobs(commands, 
                               path+'/'+data_dir+'/'+tmp_dir+'/',
                               'make_counts_files',
                               use_cluster=use_cluster)
