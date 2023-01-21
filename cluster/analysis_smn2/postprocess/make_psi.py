from post_process_input import *

if not os.path.isdir(psi_dir):
    utils.clean_dir(psi_dir)

print('\n==> Process Results Datasets')
# Load all ratios into a single dataframe
for context in ['nvs', 'dmso', 'rg']:
    min_ct = 5
    psi_df = pd.DataFrame()
    tot_df = pd.DataFrame()
    for lib_num in [1,2,3]:
        for rep_num in [1,2,3]:
            results_file = f'{results_dir}/results.smn2_{context}_lib{lib_num}_rep{rep_num}.txt'
            results_df = pd.read_csv(results_file, delim_whitespace=True)
            #print('===> Processing {}'.format(results_file))
            # Group results file by ss
            ss_df = results_df.groupby('ss').sum()
            # Filter out ss with tot < min_ct
            ix = ss_df['tot_ct'] > min_ct
            ss_df = ss_df[ix]
            # Compute ratios
            ss_df['ratio'] = ss_df['inc_ct']/ss_df['tot_ct']
            cons1_df = results_df[results_df['ss']=='ACAGGTAAGT']
            cons2_df = results_df[results_df['ss']=='CCAGGTAAGT']
            cons3_df = results_df[results_df['ss']=='GCAGGTAAGT']
            cons4_df = results_df[results_df['ss']=='TCAGGTAAGT']
            cons1_ratio = cons1_df['inc_ct'].sum() / cons1_df['tot_ct'].sum()
            cons2_ratio = cons2_df['inc_ct'].sum() / cons2_df['tot_ct'].sum()
            cons3_ratio = cons3_df['inc_ct'].sum() / cons3_df['tot_ct'].sum()
            cons4_ratio = cons4_df['inc_ct'].sum() / cons4_df['tot_ct'].sum()
            # Merge into psi_df
            cons_ratio = 0.25*(cons1_ratio+cons2_ratio+cons3_ratio+cons4_ratio)
            col = f'smn2_select_lib{lib_num}_rep{rep_num}'
            psi_df[col] = (100/cons_ratio)*ss_df['ratio'].copy()
            col = f'smn2_select_lib{lib_num}_rep{rep_num}'
            tot_df[col] = ss_df['tot_ct']
    # Filter the splice sites which are not in the order sites
    splice_df = pd.read_csv(dir_path+'/ordered_splice_sites.txt', delim_whitespace=True, header=None)
    splice_df.columns=['splice_site']
    psi_df = psi_df[psi_df.index.isin(set(splice_df['splice_site']))]
    tot_df = tot_df[tot_df.index.isin(set(splice_df['splice_site']))]
    psi_file = psi_dir+f'/psi_smn2_{context}.csv'
    tot_file = psi_dir+f'/total_smn2_{context}.csv'
    print(f'\n ==> writing the psi file {psi_file} \n')
    psi_df.to_csv(psi_file, sep=',', na_rep='NaN')
    tot_df.to_csv(tot_file, sep=',', na_rep='NaN')
