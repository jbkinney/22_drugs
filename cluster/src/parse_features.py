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
LIDs_sheet = metadata_class.excel_file['LIDs']
LIDs_sheet['USE'] = LIDs_sheet['USE'].astype(str)
LIDs_sheet = LIDs_sheet[LIDs_sheet['USE']=='True']
sample_sheet=metadata_class.excel_file['samples']
path = os.getcwd()
# Intermediate directory where split files are
interm_split_dir = path+'/'+data_dir+'/'+interm_dir+'/split_fatq'

# Get the LIDs from the LID and sample sheets
LIDs      = LIDs_sheet['LID']
LIDs_smaple_sheet = sample_sheet['LID']

# Get the index of each LID in sample sheet
# And create a dictionary of each LID with corresponding
# index in the sample sheet
# Create LID dictionary.
# LID_dict is the nested dictionary.
# In the top level the key is the LID.
# For each LID we have a dictionary of vlues: sample_name, barcode,
# read1 and read2
# If read1 or read2 is false the list of file is the keyword False
# If read1 or read2 is true we have the list of the split fastq as the 
# values of read1 and read2 keys
# read2 list is False as the same of the lenght of read1 list.

LID_dict={}
for LID in LIDs:
    idx = sample_sheet.index[LIDs_smaple_sheet==LID]
    LID_dict[LID] = {}
    LID_dict[LID]['barcode'] = sample_sheet['barcode'][idx]
    LID_dict[LID]['sample_name']  = sample_sheet['sample'][idx]
    for i in idx:
        use_read1 = sample_sheet['use_read1'][i]
        use_read2 = sample_sheet['use_read2'][i]
        if use_read1==True:
            out_file_prefix= interm_split_dir+'/read1.'+str(LID)+'.fastq.*'
            LID_dict[LID]['read1']=sorted(glob.glob(out_file_prefix))
            m = len(LID_dict[LID]['read1'])
            if (m == 0):
                raise Exception('Read1 is true but there is no split_fastq file for that')
        else:
            LID_dict[LID]['read1']=False

        if use_read2==True:
            out_file_prefix= interm_split_dir+'/read2.'+str(LID)+'.fastq.*'
            LID_dict[LID]['read2']=sorted(glob.glob(out_file_prefix))
            m = len(LID_dict[LID]['read2'])
            if (m == 0):
                raise Exception('Read1 is true but there is no split_fastq file for that')
        else:
            LID_dict[LID]['read2']=[False]

        if (use_read1==True and use_read2!=True):
            LID_dict[LID]['read2']=[False]*m


commands = []        
for LID in LIDs:
    i = 0
    for read1 in LID_dict[LID]['read1']:
        read2 = LID_dict[LID]['read2'][i]
        # The following is a python command to
        # make feature files.
        #subprocess.call(['python',
        #                 path+'/src/make_features_files.py',
        #                 '-d', path+'/'+data_dir,
        #                 '-l', str(LID),
        #                 '-r1', str(read1),
        #                 '-r2', str(read2)])

        command = 'python %s -d %s -l %s -r1 %s -r2 %s'% \
                  (path+'/src/make_features_files.py',
                   path+'/'+data_dir,
                   str(LID),
                   str(read1),
                   str(read2))
        commands.append(command)
        i = i+1

# Submit commands
utils.submit_and_complete_jobs(commands,
                               path+'/'+data_dir+'/'+tmp_dir+'/',
                               'make_feature_files',
                               use_cluster=use_cluster)
