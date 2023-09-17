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
            #if verbose:
            #    print(f'In loop {num_loops}')

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
                            #if verbose:
                            #    print(f'Motif altered to add {b} at {pos:+d}')
                                
        #if verbose:
        #    print(f'Final motif found in {num_loops} loops: {motif}')
        
        # Save motif
        motif_df.loc[trial_num,'motif'] = motif

    # Collapse motif_df so motifs are unique, with counts listed
    motif_df['trials'] = 1
    if verbose:
        print(motif_df.value_counts())
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

# Returns coordinates of positive and negative set polygons
def get_cut_polygons(pos_fold, 
                     neg_fold,
                     lo_psi, 
                     hi_psi,
                     bg_x_psi,
                     bg_y_psi,
                     lims):
    
    # Create curved part for positive selection
    c_grid = np.logspace(-4, 4, 1000)
    
    r_x = c_grid + bg_x_psi/100
    r_pos_y = pos_fold*c_grid + bg_y_psi/100
    r_neg_y = neg_fold*c_grid + bg_y_psi/100
    r_inv_y = (1/neg_fold)*c_grid + bg_y_psi/100
    
    psi_x = 100*r_x/(1+r_x)
    psi_pos_y = 100*r_pos_y/(1+r_pos_y)
    psi_neg_y = 100*r_neg_y/(1+r_neg_y)
    psi_inv_y = 100*r_inv_y/(1+r_inv_y)
    
    ix = (psi_pos_y > lo_psi) & (psi_x < hi_psi)
    pos_curve_xy = zip(psi_x[ix], psi_pos_y[ix])
    
    ix = (psi_neg_y > lo_psi) & (psi_x < hi_psi)
    neg_curve_xy = zip(psi_x[ix], psi_neg_y[ix])
    
    ix = (psi_x > lo_psi) & (psi_inv_y < hi_psi)
    inv_curve_xy = zip(psi_x[ix], psi_inv_y[ix])

    xmin = lims[0]
    xmax = lims[1]
    ymin = lims[0]
    ymax = lims[1]
    
    pos_xy = [(xmin, lo_psi)] + \
            list(pos_curve_xy) + [
            (hi_psi, ymax),
            (xmin, ymax),
            (xmin, lo_psi)]
    
    neg_xy = [(lo_psi, lo_psi)] + \
             list(neg_curve_xy) + \
             [(hi_psi, hi_psi)] +\
             list(inv_curve_xy)[::-1] +\
             [(lo_psi, lo_psi)]
    
    return pos_xy, neg_xy

# Returns array of bools stating whether polygon defined by points polygon_xy contains the points in points_xy
def does_polygon_contain(polygon_xy, points_xy):
    from shapely.geometry import Point, Polygon
    poly = Polygon(polygon_xy)
    points = [Point([x,y]) for x,y in points_xy]
    return np.array([poly.contains(point) for point in points])