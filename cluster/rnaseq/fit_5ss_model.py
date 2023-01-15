import sys
from csv import DictWriter
import pandas as pd
from AS_quant.models.regulation import ManifoldModel


def get_exon_id(event_id):
    chrom, coords = event_id.split(':')
    strand, coords = coords[-1], coords[:-1]
    exons = coords.split(',')
    return('{}:{}{}'.format(chrom, exons[1], strand))


def get_exon_intron_id(event_id, intron):
    chrom, coords = event_id.split(':')
    strand, coords = coords[-1], coords[:-1]
    exons = [x.split('-') for x in coords.split(',')]
    introns = [(e1, s2) for (_, e1), (s2, _) in zip(exons, exons[1:])]
    return('{}:{}{}'.format(chrom, introns[intron], strand))


def get_independent_exons(df):
    ids = pd.DataFrame({'exon': [get_exon_id(event_id) for event_id in df.index],
                        'i1': [get_exon_intron_id(event_id, 0) for event_id in df.index],
                        'i2': [get_exon_intron_id(event_id, 1) for event_id in df.index],
                        'event': df.index})
    ids.drop_duplicates('exon', inplace=True)
    ids.drop_duplicates('i1', inplace=True)
    ids.drop_duplicates('i2', inplace=True)
    event_ids = ids['event'].values
    return(event_ids)


def select_ids(dfs, ids):
    new_dfs = [df.loc[ids, :] for df in dfs]
    return(new_dfs)


def load_counts(prefix):
    fpath = '{}.inclusion.csv'.format(prefix)
    inclusion = pd.read_csv(fpath, index_col=0)
    
    fpath = '{}.total.csv'.format(prefix)
    total = pd.read_csv(fpath, index_col=0)
    
    fpath = '{}.log_bias.csv'.format(prefix)
    log_bias = pd.read_csv(fpath, index_col=0)
    
    ids = get_independent_exons(total)
    inclusion, total, log_bias = select_ids([inclusion, total, log_bias], ids)
    
    return(inclusion, total, log_bias)


if __name__ == '__main__':
    min_exons = 10
    conditions = ['branaplam', 'risdiplam']

    prefix = sys.argv[1]
    seqs_5ss_fpath = sys.argv[2]
    design_fpath = sys.argv[3]
    out_fpath = sys.argv[4]
    
    design = pd.read_csv(design_fpath, index_col=0)
    counts = load_counts(prefix)
    seqs_5ss = pd.read_csv(seqs_5ss_fpath).set_index('event_id')
    seqs_5ss = seqs_5ss.loc[counts[0].index, :].reset_index()
    seqs_5ss.columns = ['event_id', '5ss']
    print(seqs_5ss)
    
    with open(out_fpath, 'w') as fhand:
        fieldnames = ['seq', 'drug', 'E[beta]', 'beta_2.5', 'beta_97.5']
        writer = DictWriter(fhand, fieldnames=fieldnames)
        writer.writeheader()
        
        for seq, event_ids in seqs_5ss.groupby(['5ss'])['event_id']:
            if event_ids.shape[0] < min_exons:
                continue
            inclusion, total, log_bias = select_ids(counts, event_ids)
            model = ManifoldModel(design, recompile=False)
            model.observations.load(inclusion, total, log_bias)
            model.fit()
            s = model.get_summary()
            s['seq'] = seq
            s['drug'] = s.index
            for record in s.to_dict(orient='index').values():
                writer.writerow(record)
            fhand.flush()
