from post_process_input import *

# If both ss and bc exists it will look like old justin cipher file.
# If just bc exists it will look like counts file with extra junk columns.

# Load metadata file
metadata_class = metadata(metadata_file)
metadata_class.load_metadata()
sample_sheet  = metadata_class.excel_file['samples']
read_sheet    = metadata_class.excel_file['reads']

if not os.path.isdir(cipher_dir):
    utils.clean_dir(cipher_dir)

for sample in sample_sheet['sample']:
    print('Make Cipher file for {}'.format(str(sample))) 
    infile = count_file_dir+'/counts.'+str(sample)+'.txt'
    df = pd.read_csv(infile, sep='\t')
    df.sort_values(by=['bc','ct'], 
            axis=0, ascending=[True, False], 
            inplace=True)
    # Fill NaN with 0
    df.fillna(0, inplace=True)
    # Group by bc. For a given bc, most prevalent ss will be the first one
    g = df.groupby('bc')
    out_df = g.first()
    # Record the total number of ss observed for each bc
    tmp_df = g.size()
    out_df.loc[tmp_df.index, 'numss'] = tmp_df
    # Record otherct, i.e. the total number of reads linking each bc to a DIFFERENT ss
    tmp_df = g.agg(np.sum)
    out_df.loc[tmp_df.index,'totct'] = tmp_df['ct']
    out_df['otherct'] = out_df['totct']-out_df['ct']
    # Compute indices that are ok to use
    # Criterion for use: ct >= 2 and ct >= 4*otherct
    ok_indices = (out_df['ct'] >= 2) & (out_df['ct'] >= 4*out_df['otherct'])
    # Get final dataframe of associated indices
    out_df = out_df[ok_indices]
    #out_df = out_df.reindex[ok_indices, ['ss','ct','otherct','numss']]
    out_df = out_df.reindex(columns=['ss','ct','otherct','numss'])
    out_df.sort_values(by='ct', axis=0, ascending=False, inplace=True)
    # Drop NAN from the file without ss
    out_df = out_df.dropna(axis=1, how='all')
    # Reset index
    out_df = out_df.reset_index()
    # Write output file
    out_file = cipher_dir+'/cipher.'+str(sample)+'.txt'
    out_df.to_csv(out_file, sep='\t')
