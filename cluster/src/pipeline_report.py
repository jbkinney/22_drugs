import utils, sys, os
import pandas as pd

def report_features(df, LID, num_reads, dir, file_extension):
    """
    Generate report on fly.
    """
    # load param file
    # Make report directory under the intermediate directory
    report_dir = dir+'/reports'
    if not os.path.isdir(report_dir):
        utils.clean_dir(report_dir)

    # Find unique sample names
    unique_sample = df['sample'].unique()
    rep_data = {}
    rep_data.setdefault('sample', [])
    rep_data.setdefault('successful reads', [])
    rep_data.setdefault('LID', LID)
    rep_data.setdefault('total number of lines', num_reads)
    
    for u in unique_sample:
        rep_data['sample'].append(u)
        rep_data['successful reads'].append((df['sample'].values==u).sum())
        
    report_df = pd.DataFrame(rep_data, 
                             columns=['sample', 
                                      'LID', 
                                      'successful reads',
                                      'total number of lines'])
    out_file=str(report_dir+'/report.'+str(LID)+file_extension+'.txt')
    report_df.to_csv(out_file, sep='\t', index=False, na_rep='None')
