import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import numpy as np
import pandas as pd

configs = [{'background': 'smn2_pt1', 
            'drug1': 'Risdiplam', 
            'drug2': 'ASOi6', 
            'sheet_names': ['smn2_pt1_ris', 
                            'smn2_pt1_asoi6', 
                            'smn2_pt1_ris_asoi6']},
           {'background': 'smn2_pt1', 
            'drug1': 'Risdiplam', 
            'drug2': 'ASOi7', 
            'sheet_names': ['smn2_pt1_ris', 
                            'smn2_pt1_asoi7_v3', 
                            'smn2_pt1_ris_asoi7_v2']},
           {'background': 'smn2_pt1', 
            'drug1': 'Risdiplam', 
            'drug2': 'Branaplam',
            'sheet_names': ['smn2_pt1_ris', 
                            'smn2_pt1_bran_v2', 
                            'smn2_pt1_ris_bran']},
           {'background': 'elp1', 
            'drug1': 'Rectas', 
            'drug2': 'ASO',
            'sheet_names': ['elp1_rectas', 
                            'elp1_aso', 
                            'elp1_rectas_aso']},
           {'background': 't25a', 
            'drug1': 'Risdiplam', 
            'drug2': 'ASOi7',
            'sheet_names': ['t25a_ris', 
                            't25a_asoi7', 
                            't25a_ris_asoi7']},
           {'background': 't25a', 
            'drug1': 'Risdiplam', 
            'drug2': 'ASOi6',
            'sheet_names': ['t25a_ris', 
                            't25a_asoi6', 
                            't25a_ris_asoi6']},
           {'background': 't25a', 
            'drug1': 'Risdiplam', 
            'drug2': 'Branaplam', 
            'sheet_names': ['t25a_ris', 
                            't25a_bran', 
                            't25a_ris_bran']}]

EC2x = {'Risdiplam': 14,
        'Branaplam': 7,
        'ASOi6': 0.6,
        'ASOi7': 0.1,
        'Rectas': 0.2,
        'ASO': 0.08}


def load_data(sheet_name):
    file_name = '../data/22.08.05_qPCR_titrations.xlsx'
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

import pandas as pd
import numpy as np
import re 
import mplcursors

# Function to load data
def get_data_from_cvs_files(filenames):
    """
    Read the csv file and return the pandas dataframe with
    the median PSI values.
    """
    
    for i, filename in enumerate(filenames):
        col = filename.split('/')[-1].split('.')[0][4:]
        tmp_df = pd.read_csv(filename, sep=',')
        if 'brca2' in col:
            tmp_df['ss'] = [('A'+seq).replace('T','U') for seq in tmp_df['ss']]
        elif 'ikbkap' in col:
            tmp_df['ss'] = [('A'+seq).replace('T','U') for seq in tmp_df['ss']]
        elif 'smn2' in col:
            tmp_df['ss'] = [seq.replace('T','U') for seq in tmp_df['ss']]
        else:
            assert False, 'This should not happen'
        #print(tmp_df.columns)
        tmp_df = tmp_df.set_index('ss')
        tmp_df[col]= tmp_df.median(axis=1)
        tmp_df = tmp_df[[col]]
        if i==0:
            df = tmp_df
        else:
            df = pd.merge(left=df, right=tmp_df, left_index=True, right_index=True, how='outer')
    return df

# Define IUPAC/motif handling functions
_iupac_to_bs_dict = {
        'A':('A',),
        'C':('C',),
        'G':('G',),
        'U':('U',),
        'R':('A','G',),
        'Y':('C','U',),
        'S':('C','G',),
        'W':('A','U',),
        'K':('G','U',),
        'M':('A','C',),
        'B':('C','G','U',),
        'D':('A','G','U',),
        'H':('A','C','U',),
        'V':('A','C','G',),
        'N':('A','C','G','U',)
    }
_iupac_to_regex_dict = {k:f'[{"".join(v)}]' for k,v in _iupac_to_bs_dict.items()}
_bs_to_iupac_dict = {v: k for k, v in _iupac_to_bs_dict.items()}

# Funciton to make plots interactive
def add_click_labels(x, y, seqs, ax, ix=None):
    N = len(x)
    assert len(y)==N
    if ix is None:
        ix = np.ones(N).astype(bool)
    p1ot = ax.scatter(x[ix],y[ix], facecolor='none', edgecolor='none')
    mplcursors.cursor([p1ot],hover=False).connect(
        "add", lambda sel: sel.annotation.set_text(seqs[ix][sel.target.index]))

# Get list of nucleotides given iupac symbol
def iupac_nt_to_list(nt):
    assert(nt in _iupac_to_bs_dict.keys())
    return list(_iupac_to_bs_dict[nt])

# Function to compute regular expression from motif
def iupac_to_regex(iupac, rna=True):
    '''Given an IUPAC motif, returns an equivalent regular expresion'''
    s = iupac
    for k, v in _iupac_to_regex_dict.items():
        s = s.replace(k,v)
    if not rna:
        s = s.replace('U','T')
    return s

# Function to compute the "mass" of a motif
def motif_to_mass(motif):
    '''Computes the mass of an IUPAC motif, defined as the number of sequences that motif hits.'''
    mass = 1
    poss = [-4,-3,-2,-1,3,4,5,6] 
    for pos in poss:
        i = pos+4   # Is right if motif contains '/'
        iupac = motif[i]
        bs = _iupac_to_bs_dict[iupac]
        num_bs = len(bs)
        mass *= num_bs
    return mass

# Function to convert motif to indices
def motif_to_ix(motif, seqs):
    '''Given an IUPAC motif and a set of sequences, returns the indices of matching sequences'''
    seqs = np.array(seqs)
    regex = iupac_to_regex(motif)
    ix    = np.array([bool(re.match(regex, seq)) for seq in seqs])
    return ix

# Compute necessary and sufficient motifs for a given set of included
# and excluded sequences
def get_sufficient_and_necessary_motifs(seqs, in_ix, ex_ix, num_trials=100, verbose=False):

    # First compute most restrictive motif that matches all 'in' sequences
    # If this motif also excludes all excluded sequences, it will be the 
    # "sufficient motif"
    motif = 'xxxx/GUxxxx'
    poss = [-4,-3,-2,-1,3,4,5,6] 
    motif_list = list(motif)
    bases = list('ACGU')
    seqs = np.array(seqs)
    if in_ix is not None:
        in_seqs = seqs[in_ix]
    for pos in poss:
        i = pos+4
        bs = list(set([seq[i] for seq in in_seqs]))
        bs.sort()
        bs = tuple(bs)
        iupac = _bs_to_iupac_dict[bs]
        motif_list[i] = iupac
    all_in_motif = ''.join(motif_list)
    
    # If this motif does NOT match any excluded sequences, it is a sufficient motif
    ix = motif_to_ix(motif=all_in_motif,
                           seqs=seqs)
    if sum(ix & ex_ix)>0:
        print('Could not find sufficient motif. Returning all_in_motif instead')
        return all_in_motif, None
    
    # Define sufficient motif
    sufficient_motif = all_in_motif
                
    # Create dataframe to hold maximal motifs from multiple trials
    motif_df = pd.DataFrame(index=range(num_trials), columns=['motif'])
    
    # Perform trials
    for trial_num in range(num_trials):

        # Initialize using from sufficient motif
        motif = sufficient_motif
        motif_altered = True
        num_loops = 0
        
        # While motif keeps being altered, iterate over all pos and b
        while motif_altered:
            motif_altered = False
            num_loops += 1
            if verbose:
                print(f'In loop {num_loops}')

            # Iterate over randomized positions pos
            np.random.shuffle(poss)
            for pos in poss:
                i = pos+4

                # Iterate over randomized bases b
                np.random.shuffle(bases)
                for b in bases:

                    # If b is not already permitted by motif
                    iupac = motif[i]
                    allowed_bs = _iupac_to_bs_dict[iupac]
                    if b not in allowed_bs:

                        # Add b to allowed bs
                        new_bs = list(allowed_bs)
                        new_bs.append(b)
                        new_bs.sort()
                        new_bs = tuple(new_bs)

                        # Get iupac symbol for bs
                        new_iupac = _bs_to_iupac_dict[new_bs]

                        # Create new motif
                        new_motif_list = list(motif)
                        new_motif_list[i] = new_iupac
                        new_motif = ''.join(new_motif_list)

                        # Determine new sequences hit
                        new_re = iupac_to_regex(new_motif)
                        new_ix = np.array([bool(re.match(new_re, seq)) for seq in seqs])

                        # Compute number of non-activted seqs
                        nact_seqs_hit = sum(new_ix & ex_ix)
                        act_seqs_hit = sum(new_ix & in_ix)
                        N_act = sum(in_ix)

                        # If new motif does not hit non-activated seqs AND hits all activated sequences, save
                        if (nact_seqs_hit == 0) and (act_seqs_hit == N_act):
                            motif = new_motif
                            motif_altered = True
                            if verbose:
                                print(f'Motif altered to add {b} at {pos:+d}')
                                
        if verbose:
            print(f'Final motif found in {num_loops} loops: {motif}')
        
        # Save motif
        motif_df.loc[trial_num,'motif'] = motif

    # Collapse motif_df so motifs are unique, with counts listed
    motif_df['trials'] = 1
    if verbose:
        print(motif_df)
    motif_df = motif_df.groupby('motif').sum().reset_index()
    
    # Compute mass and num positives for each identified motif
    motif_df['mass'] = 0
    for i in range(len(motif_df)):
        motif = motif_df.loc[i,'motif']
        motif_df.loc[i,'mass'] = motif_to_mass(motif)
        
    # Sort by mass, then by ct
    motif_df.sort_values(by=['mass', 'trials'], 
                         ascending=False, 
                         inplace=True)
    motif_df.reset_index(inplace=True, drop=True)
    
    # Define necessary motif
    necessary_motif = motif_df.loc[0,'motif']
            
    # Return sufficient and necessary motifs
    return sufficient_motif, necessary_motif

