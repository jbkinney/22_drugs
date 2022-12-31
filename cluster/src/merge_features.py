import sys, os
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
path = os.getcwd()

interm_combined_dir = data_dir+'/'+interm_dir+'/combined_feature'
if not os.path.isdir(interm_combined_dir):
    utils.clean_dir(interm_combined_dir)

for LID in LIDs_sheet['LID']:
    skip = False
    feature_files_path = str(data_dir+'/'+interm_dir+
                             '/feature/feature.'+str(LID)+'*.txt')
    feature_files_list = glob.glob(feature_files_path)
    if len(feature_files_list)==0:
        print('For LID = {} feature files are not exist.'.format(LID))
        skip = True
    if skip==False:
        result = str(interm_combined_dir+'/features_combined_'+str(LID)+'.txt')
        with open(result, 'wb') as out_file:
            for f in feature_files_list:
                with open(f, 'rb') as in_file:
                    out_file.write(in_file.read())
