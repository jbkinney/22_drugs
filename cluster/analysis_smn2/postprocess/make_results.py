from post_process_input import *

# Load metadata file
metadata_class = metadata(metadata_file)
metadata_class.load_metadata()
LIDs_sheet  = metadata_class.excel_file['LIDs']
LIDs_sheet['USE'] = LIDs_sheet['USE'].astype(str)
LIDs_sheet   = LIDs_sheet[LIDs_sheet['USE']=='True']
list_lids    = LIDs_sheet['LID'].tolist()
sample_sheet = metadata_class.excel_file['samples']
# Only USE=True LID from sample sheet is reading
sample_sheet = sample_sheet[sample_sheet['LID'].isin(list_lids)]
read_sheet   = metadata_class.excel_file['reads']

post_process = metadata_class.excel_file['postprocess']

if not os.path.isdir(results_dir):
    utils.clean_dir(results_dir)

experiment  = post_process['experiment']
sample_type = post_process['sample_type']

for unique_exp in experiment.unique():
    files = post_process['file'][experiment==unique_exp]
    # Find the library file
    library_file= files[sample_type=='cipher'].values[0]
    library_file=cipher_dir+'/'+library_file
    library_df  = pd.read_csv(library_file, delim_whitespace=True, header=0)
    library_df  = library_df[['bc','ss','ct','otherct']].\
            rename(columns={'ct':'lib_ct','otherct':'mis_ct'})

    print('--> Processing the experiment {}'.format(unique_exp)) 
    # Find the barcode libraries
    bc_files     = cipher_dir+'/'+files[sample_type!='cipher']
    bc_types     = sample_type[bc_files.index]
    out_string =  'library_file:\t%s\n'%unique_exp
    bc_dfs = []
    cols = ['lib_ct', 'mis_ct']
   
    for _bc in bc_files:
        bc_name =  str(bc_types[bc_files==_bc].values[0])
        ct_col_label = str(bc_name+ '_ct')
        bc_df = pd.read_csv(_bc, delim_whitespace=True, header=0)
        cols.append(ct_col_label)
        bc_df = bc_df[['bc','ct']].\
                rename(columns={'ct':ct_col_label})
        bc_dfs.append(bc_df)
        out_string += '%s_file:\t%s\n'%(bc_name, _bc)
    
    results_df = library_df.copy()

    for df in bc_dfs:
        results_df = results_df.merge(df, how='left', on='bc')
    
    results_df.fillna(0,inplace=True)
    results_df.sort_values(by='lib_ct',inplace=True,ascending=False)
    results_df.reset_index(inplace=True,drop=True)
    
    out_file = results_dir+'/results.'+str(unique_exp)+'.txt'
    results_df.to_csv(out_file, sep='\t', float_format='%d', na_rep='NaN')
