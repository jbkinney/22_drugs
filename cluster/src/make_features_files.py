# This file will run on computational nodes.
import utils,sys,os
from metadata import *
import argparse
from pipeline_report import *

parser = argparse.ArgumentParser(description='Making feature files')

parser.add_argument('-d','--data_dir', 
                    help='Project metadata directory', 
                    required=True)

parser.add_argument('-l','--LID',
                    help='LID value as string', 
                    required=True)

parser.add_argument('-r1','--read1',
                    help='Splitted read1 file, abs path included', 
                    required=True)

parser.add_argument('-r2','--read2',
                    help='Splitted read2 file, abs path included', 
                    required=True)

args = parser.parse_args()

data_dir = args.data_dir
LID = int(args.LID)
read1_file = args.read1
read2_file = args.read2

# load metadata file
sys.path.append(data_dir)
from params import *
metadata_class = metadata(data_dir+'/'+metadata_file)
metadata_class.load_metadata()
sample_sheet  = metadata_class.excel_file['samples']
read_sheet    = metadata_class.excel_file['reads']
LIDs_in_sample_sheet = sample_sheet['LID']
# Get extension name of the read1 file
filename, file_extension = os.path.splitext(read1_file)

# Make feature directory under the intermediate directory
interm_feature_dir = data_dir+'/'+interm_dir+'/feature'
if not os.path.isdir(interm_feature_dir):
    utils.clean_dir(interm_feature_dir)

# Find all rows from sample sheet which has 
# current LID on computational node
idx = sample_sheet.index[LIDs_in_sample_sheet==LID]
barcodes            = sample_sheet['barcode'][idx]
sample_name         = sample_sheet['sample'][idx]
# tolist() will change to to_list() in feature release of pandas

# Make search dictionary search_dict = {barcode:sample_name}
search_dict = dict(zip(barcodes, sample_name))

# Open read files
r1_f = open(read1_file)
# If read2 file is not passed (false), r2_f = False
if read2_file!='False':
    r2_f = open(read2_file)
else:
    r2_f = False

# Make sample_dict contains regx
sample_dict={}
for sample in sample_name:

    idx = sample_sheet.index[(sample_sheet['LID']==LID) & (sample_sheet['sample']==sample)].tolist()
    # Get the name of the amplicon for read1 of current sample
    read1_amp = sample_sheet['read1'][idx]
    # Find the regex index corresponding to read1_amp
    jdx = list(np.where(read_sheet['read'].values==read1_amp.values)[0])

    sample_dict[sample]={}
    sample_dict[sample]['read1']=read_sheet['regex'][jdx].tolist()[0]
    
    if r2_f:
        read2_amp = sample_sheet['read2'][idx]
        jdx = list(np.where(read_sheet['read'].values==read2_amp.values)[0])
        sample_dict[sample]['read2']=read_sheet['regex'][jdx].tolist()[0]

    # Take union of the features
    features_col       = sample_sheet['features'][idx].tolist()
    features_item = []
    for item in features_col:
        features_item.append(item.split(','))
    features_list = list(set().union(*features_item))
    # Remove whitespaces
    features_list = [x.strip(' ') for x in features_list]
    
    sample_dict[sample]['features_list']=features_list



def match_barcode(seq, search_dict, search_area=20):
    tag_length = 0
    region = False
    for barcode in search_dict.keys():
        k = seq.find(barcode, 0, search_area)
        if k >= 0:
            region = search_dict[barcode]
            tag_length = len(barcode)+k
    return (region, tag_length)


def parse_read(sequence, regex, feature_list):
    m = re.match(regex,sequence)
    # If parse failed, return False
    if not m:
        m = False
        write_result = False
        out_dict = {}
        for key in features_list:
            out_dict[key] = None
    # Otherwise, return out_dict with all possible elements
    else:
        raw_dict = m.groupdict()
        out_dict = {}
        for key in features_list:
            key_rc = key+'_rc'
            if key in raw_dict:
                out_dict[key]=str(raw_dict[key])
                write_result  = True
            elif key_rc in raw_dict:
                out_dict[key]=utils.rc(str(raw_dict[key_rc]))
                write_result = True
            else:
                out_dict[key]= None
                write_result = False
    return m, out_dict, write_result


stop = False
num_reads = 0
# Main loop over all lines of the splitted read files
# st_results: Structured list of results
st_results = []

while not stop:

    # Get reads; halt loop if reads have length 0
    read1 = utils.get_next_read_from_fastq(r1_f)
    # If reached end of the read1 file stop
    if len(read1) == 0:
        stop = True
        continue

    # Adopt to the condition when read2 exists
    if r2_f:
        read2 = utils.get_next_read_from_fastq(r2_f)
        if len(read2)==0:
            stop=True
            continue
    num_reads += 1

    # Determine sequence sample by matching barcode, and clip barcode
    sample, tag_length = match_barcode(read1,
                                       search_dict,
                                       search_area=20)
    assert len(read1) > tag_length, 'Read 1 is not long enough'

    # Clip barcode from read2 if it exists
    if r2_f:
        sample2, tag_length2 = match_barcode(read2,
                                             search_dict,
                                             search_area=20)
        assert len(read2) > tag_length2, 'Read 2 is not long enough'
    else:
        sample2 = False

    # Barcode did not found in read1 but it is available in read2
    if sample==False:
        write_result = False
        m2 = False
        if sample2!=False:
            read2_crop=str(read2[tag_length2:])
            regx2 = re.compile(sample_dict[sample2]['read2'])
            features_list=sample_dict[sample2]['features_list']
            m2, out_dict, write_result = parse_read(read2_crop, regx2, features_list)

        if write_result:
            out_dict['sample'] = sample2
            st_results.append(out_dict) 
         
    if sample!=False:
        write_result = False
        # Crop the barcode from the sequence
        read1_crop = str(read1[tag_length:])
        regx1 = re.compile(sample_dict[sample]['read1'])
        features_list=sample_dict[sample]['features_list']
        m1, out_dict1, write_result1 = parse_read(read1_crop, regx1, features_list)
        # Check the read2 region:

        if sample2!=False:
           read2_crop=str(read2[tag_length2:])
           regx2 = re.compile(sample_dict[sample2]['read2'])
           features_list=sample_dict[sample2]['features_list']
           m2, out_dict2, write_result2 = parse_read(read2_crop, regx2, features_list)
        
        # Read2 file is not exists or there was not region2 -> Sample2 == False
        # We just work with read1 file
        if sample2==False:
            write_result2 = False
            m2 = False
            out_dict2 = {}
            for key in features_list:
                out_dict2[key] = None
             
        out_dict = {}
        for key in features_list:
            f1 = out_dict1[key]
            f2 = out_dict2[key]
            if (f1 and f2):
                if (f1!=f2): 
                    write_result = False
                elif(f1==f2 and f1!='None'):
                    write_result = write_result1
                    out_dict[key] = f1
            elif f1:
                write_result = write_result1
                out_dict[key] = f1
            elif f2:
                write_result = write_result2
                out_dict[key] = f2
        # append all successful results in structured list
        if write_result:
            out_dict['sample'] = sample
            st_results.append(out_dict)
# The output has #features+1 numbers of columns
out_file=str(interm_feature_dir+'/feature.'+str(LID)+file_extension+'.txt')

cols = ['sample']
for key in features_list:
    cols.append(key)
df_res = pd.DataFrame(st_results, columns=cols)

# Check for all features REQUESTED we dont have any NONE
# This is important where we have ss and bc in read1 and read2 separately.
for key in features_list:
    df_res = df_res[df_res[key].notna()]

df_res.to_csv(out_file, sep='\t', index=False, na_rep='None')

# Make report
report_features(df_res, LID, num_reads, data_dir+'/'+interm_dir, file_extension)
