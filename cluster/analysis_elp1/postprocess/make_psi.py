from post_process_input import *

# Names in the manuscript is different
# ikbkap -> elp

if not os.path.isdir(psi_dir):
    utils.clean_dir(psi_dir)

drug_name = ['dmso', 'rg', 'nvs']

print('\n==> Process Results Datasets')
# Load all ratios into a single dataframe
for k, context in enumerate(['dmso', 'rg', 'nvssm']):
    min_ct = 2
    psi_df = pd.DataFrame()
    tot_df = pd.DataFrame()
    for lib_num in [1,2]:
        for rep_num in [1,2]:
            results_file = f'{results_dir}/results.ikbkap_{context}_lib{lib_num}_rep{rep_num}.txt'
            results_df = pd.read_csv(results_file, delim_whitespace=True)
            #print('===> Processing {}'.format(results_file))
            # Group results file by ss
            ss_df = results_df.groupby('ss').sum()
            # Filter out ss with tot < min_ct
            ix = ss_df['tot_ct'] > min_ct
            ss_df = ss_df[ix]
            # Compute ratios
            ss_df['ratio'] = ss_df['inc_ct']/ss_df['tot_ct']
            cons_df = results_df[results_df['ss']=='CAAGTAAGT']
            cons_ratio = cons_df['inc_ct'].sum() / cons_df['tot_ct'].sum()
            col = f'ikbkap_select_lib{lib_num}_rep{rep_num}'
            psi_df[col] = (100/cons_ratio)*ss_df['ratio'].copy()
            col = f'ikbkap_select_lib{lib_num}_rep{rep_num}'
            tot_df[col] = ss_df['tot_ct']
    drug = drug_name[k]
    psi_file = psi_dir+f'/psi_elp1_{drug}.csv'
    tot_file = psi_dir+f'/total_elp1_{drug}.csv'
    print(f'\n ==> writing the psi file {psi_file} \n')
    psi_df.to_csv(psi_file, sep=',', na_rep='NaN')
    tot_df.to_csv(tot_file, sep=',', na_rep='NaN')
