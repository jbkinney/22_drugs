# This file will run on computational nodes.
import utils,sys,os
from metadata import *
import argparse

parser = argparse.ArgumentParser(description='Making counts files')

parser.add_argument('-d','--data_dir', 
                    help='Project metadata directory', 
                    required=True)

parser.add_argument('-s','--sample_name',
                    help='Sample name as string', 
                    required=True)

args = parser.parse_args()
data_dir = args.data_dir
sample = args.sample_name

# load metadata file
sys.path.append(data_dir)
from params import *
metadata_class = metadata(data_dir+'/'+metadata_file)
metadata_class.load_metadata()
LIDs_sheet  = metadata_class.excel_file['LIDs']
LIDs_sheet['USE'] = LIDs_sheet['USE'].astype(str)
LIDs_sheet   = LIDs_sheet[LIDs_sheet['USE']=='True']
list_lids    = LIDs_sheet['LID'].tolist()
sample_sheet = metadata_class.excel_file['samples']
# Only USE=True LID from sample sheet is reading
sample_sheet = sample_sheet[sample_sheet['LID'].isin(list_lids)]


# Get LID corresponding to the sample_name for specific LID
# passed to this file from sample_sheet
idx = sample_sheet.index[sample_sheet['sample']==sample]
# MahdiK same sample name for different LIDs
LID = sample_sheet['LID'][idx].values[0]

# Intermediate directory where feature files are
path = os.getcwd()
interm_feature_dir = data_dir+'/'+interm_dir+'/feature'
# Make counts directory under the intermediate directory
interm_counts_dir = data_dir+'/'+interm_dir+'/counts'
if not os.path.isdir(interm_counts_dir):
    utils.clean_dir(interm_counts_dir)


# loop over feature files with LID correspond to sample name.
in_file_prefix= interm_feature_dir+'/feature.'+str(LID)+'.*'+'.txt'
in_file_list = glob.glob(in_file_prefix)

tmp_df = pd.DataFrame()
for in_file in in_file_list:
    # Load input file as data frame
    in_df = pd.read_csv(in_file, sep='\t')
    tmp_df = tmp_df.append(in_df[in_df['sample']==sample],
                           ignore_index=True,
                           sort=False)

# Get the column name, ignore the 'sample' name
cols = list(tmp_df)
cols.remove('sample')
# Count the exact occurence in rows
out_df = tmp_df.groupby(cols).size().to_frame(name = 'ct').reset_index()
out_df.sort_values(by='ct',inplace=True,ascending=False)
out_df.reset_index(inplace=True,drop=True)
# Reorder columns to bring ct in the first column
cols = out_df.columns.tolist()[::-1]
out_df  = out_df[cols]
# Write out count files
out_file=str(interm_counts_dir+'/counts.'+str(sample)+'.txt')
out_df.to_csv(out_file, sep = '\t', index=False)
