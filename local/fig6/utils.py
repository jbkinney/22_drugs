import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import numpy as np
import pandas as pd

configs = [{'background': 'smn2_pt1', 'drug1': 'Risdiplam', 'drug2': 'ASOi6', 'sheet_names': ['smn2_pt1_ris', 'smn2_pt1_asoi6', 'smn2_pt1_ris_asoi6']},
           {'background': 'smn2_pt1', 'drug1': 'Risdiplam', 'drug2': 'ASOi7', 'sheet_names': [
               'smn2_pt1_ris', 'smn2_pt1_asoi7_v3', 'smn2_pt1_ris_asoi7_v2']},
           {'background': 'smn2_pt1', 'drug1': 'Risdiplam', 'drug2': 'Branaplam',
               'sheet_names': ['smn2_pt1_ris', 'smn2_pt1_bran_v2', 'smn2_pt1_ris_bran']},
           {'background': 'elp1', 'drug1': 'Rectas', 'drug2': 'ASO',
               'sheet_names': ['elp1_rectas', 'elp1_aso', 'elp1_rectas_aso']},
           {'background': 't25a', 'drug1': 'Risdiplam', 'drug2': 'ASOi7',
               'sheet_names': ['t25a_ris', 't25a_asoi7', 't25a_ris_asoi7']},
           {'background': 't25a', 'drug1': 'Risdiplam', 'drug2': 'ASOi6',
               'sheet_names': ['t25a_ris', 't25a_asoi6', 't25a_ris_asoi6']},
           {'background': 't25a', 'drug1': 'Risdiplam', 'drug2': 'Branaplam', 'sheet_names': ['t25a_ris', 't25a_bran', 't25a_ris_bran']}]

EC2x = {'Risdiplam': 14,
        'Branaplam': 7,
        'ASOi6': 0.6,
        'ASOi7': 0.1,
        'Rectas': 0.2,
        'ASO': 0.08}


def load_data(sheet_name):
    file_name = '../analysis/22.05.20_qPCR_titration_analysis/22.05.20_qPCR_titrations.xlsx'
    data_df = pd.read_excel(file_name,
                            sheet_name=sheet_name, header=[0, 1], index_col=[0, 1])
    data_df.index.rename(['conc', 'bio_rep'], inplace=True)
    data_df.columns.rename(['primers', 'tech_rep'], inplace=True)

    # Pivot to make tidy data
    tidy_df = data_df.melt(value_name='cycles',
                           ignore_index=False).reset_index()

    df = tidy_df.pivot(index=['conc', 'bio_rep', 'tech_rep'],
                       columns='primers', values='cycles').reset_index()
    df['dCt'] = df['exclusion'] - df['inclusion']
    df['sample_id'] = ['{}.{}'.format(c, b)
                       for c, b in zip(df['conc'], df['bio_rep'])]
    df['run_id'] = ['{}.{}'.format(s, r)
                    for s, r in zip(df['sample_id'], df['tech_rep'])]
    df = df.groupby(['sample_id'])['conc', 'bio_rep', 'dCt'].mean()
    return df


def load_combination_data(config):
    drug1 = load_data(config['sheet_names'][0])
    drug2 = load_data(config['sheet_names'][1])
    mix = load_data(config['sheet_names'][2])

    drug1['conc1'] = drug1['conc'].values
    drug1['conc2'] = 0
    drug1['drug1'] = 1
    drug1['drug2'] = 0
    drug1['mix'] = 0
    drug1['drug'] = 'drug1'
    drug1['conc'] = drug1['conc'].values / EC2x[config['drug1']]

    drug2['conc1'] = 0
    drug2['conc2'] = drug2['conc'].values
    drug2['drug1'] = 0
    drug2['drug2'] = 1
    drug2['mix'] = 0
    drug2['drug'] = 'drug2'
    drug2['conc'] = drug2['conc'].values / EC2x[config['drug2']]

    mix['conc'] = mix['conc']
    mix['conc1'] = mix['conc'] * EC2x[config['drug1']] / 2
    mix['conc2'] = mix['conc'] * EC2x[config['drug2']] / 2
    mix['drug1'] = 1
    mix['drug2'] = 1
    mix['mix'] = 1
    mix['drug'] = 'mix'

    full_data = pd.concat([drug1, drug2, mix], axis=0)
    full_data['psi'] = 2**full_data['dCt'] / (1 + 2**full_data['dCt'])
    full_data['total_conc'] = full_data['conc1'] + full_data['conc2']
    return full_data


def get_numpyro_data(df, conc, drugs):
    sel_idxs = [x in drugs for x in df['drug']]
    df = df.loc[sel_idxs, :]
    data = {}
    for conc_val in conc:
        data[conc_val] = df[conc_val].values
    data['dCt'] = df['dCt'].values
    data['psi'] = df['psi'].values
    for drug_val in drugs:
        if drug_val != 'mix':
            data[drug_val] = df[drug_val].astype('bool')
    data_df = pd.DataFrame(data)
    data_df.dropna(inplace=True)
    data_df.reset_index(inplace=True)
    return data_df


def x_to_ohe(x,
             alphabet_name,
             ravel_seqs=True):
    """
    Convert a sequence array to a one-hot encoded matrix.

    Parameters
    ----------
    x: (np.ndarray)
        (N,) array of input sequences, each of length L

    alphabet_name: (np.ndarray)
        (C,) array describing the alphabet sequences are drawn from.

    ravel_seqs: (bool)
        Whether to return an (N, L*C) array, as opposed to an (N, L, C) array.

    Returns
    -------
    x_ohe: (np.ndarray)
        Array of one-hot encoded sequences, stored as np.int8 values.
    """

    # Define built-in alphabets
    alphabet_dict = {
        'dna': np.array(['A', 'C', 'G', 'T']),
        'rna': np.array(['A', 'C', 'G', 'U']),
        'protein': np.array(['A', 'C', 'D', 'E', 'F',
                             'G', 'H', 'I', 'K', 'L',
                             'M', 'N', 'P', 'Q', 'R',
                             'S', 'T', 'V', 'W', 'Y']),
        'protein*': np.array(['A', 'C', 'D', 'E', 'F',
                              'G', 'H', 'I', 'K', 'L',
                              'M', 'N', 'P', 'Q', 'R',
                              'S', 'T', 'V', 'W', 'Y', '*'])
    }

    alphabet = alphabet_dict[alphabet_name]
    # Get dimensions
    L = len(x[0])
    N = len(x)
    C = len(alphabet)

    # Shape sequences as array of int8s
    x_arr = np.frombuffer(bytes(''.join(x), 'utf-8'),
                          np.int8, N * L).reshape([N, L])

    # Create alphabet as array of int8s
    alphabet_arr = np.frombuffer(bytes(''.join(alphabet), 'utf-8'),
                                 np.int8, C)

    # Compute (N,L,C) grid of one-hot encoded values
    x_nlc = (x_arr[:, :, np.newaxis] ==
             alphabet_arr[np.newaxis, np.newaxis, :]).astype(np.int8)

    # Ravel if requested
    if ravel_seqs:
        x_ohe = x_nlc.reshape([N, L * C])
    else:
        x_ohe = x_nlc

    return x_ohe
