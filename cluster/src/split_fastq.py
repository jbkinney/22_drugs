import sys, os
from metadata import *
import utils

data_dir = sys.argv[1]
sys.path.append(data_dir)
# Import parameters from params.py in projectx directory
from params import *

# Make split_fastq directory under the intermidiate directory
interm_split_dir = data_dir+'/'+interm_dir+'/split_fatq'
if not os.path.isdir(interm_split_dir):
    utils.clean_dir(interm_split_dir)


num_lines = int(4*reads_per_split)

# Load metadata file
metadata_class = metadata(data_dir+'/'+metadata_file)
metadata_class.load_metadata()
LIDs_sheet = metadata_class.excel_file['LIDs']
path = os.getcwd()

commands = []
split_fastq_files = []
LIDs_sheet['USE'] = LIDs_sheet['USE'].astype(str)
LIDs_sheet = LIDs_sheet[LIDs_sheet['USE']=='True']
LIDs      = LIDs_sheet['LID']

for idx in LIDs.index:
    LID = LIDs[idx]
    read1 = LIDs_sheet['read1_file'][idx]
    read2 = LIDs_sheet['read2_file'][idx]
    if read1!='NONE':
        file = illumina_run_dir+'/'+read1
        out_file_prefix=str(path+'/'+interm_split_dir+'/read1.'+str(LID)+'.fastq.')
        command = 'zcat %s | split -a 3 -l %d - %s'%\
                  (file,num_lines,out_file_prefix)
        split_fastq_files.append(out_file_prefix)
        commands.append(command)
    if read2!='NONE':
        file = illumina_run_dir+'/'+read2
        out_file_prefix=str(path+'/'+interm_split_dir+'/read2.'+str(LID)+'.fastq.')
        command = 'zcat %s | split -a 3 -l %d - %s'%\
                  (file,num_lines,out_file_prefix)
        split_fastq_files.append(out_file_prefix)
        commands.append(command)

utils.submit_and_complete_jobs(commands,
                               path+'/'+data_dir+'/'+tmp_dir+'/',
                               'make_split_fastq_files',
                               use_cluster=use_cluster)
