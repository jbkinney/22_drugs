import sys
import pandas as pd

from csv import DictWriter
from os.path import join, abspath, dirname
from tqdm import tqdm

from pysam import Fastafile

from AS_quant.genome import Exon


def get_5ss_sequences(event_ids, fastafile):
    for event_id in event_ids:
        chrom, coords = event_id.split(':')
        coords, strand = coords[:-1], coords[-1]
        exons = [[int(x) for x in exon.split('-')] for exon in coords.split(',')]
        exon = Exon(chrom.lstrip('chr'), [exons[1][0]], [exons[1][1]], strand)
        seq = exon.get_5ss_seqs(fastafile)[0]
        record = {'event_id': event_id, '5ss': seq}
        yield(record)


if __name__ == '__main__':
    #genome_fpath = join(cwd, 'reference', 'Homo_sapiens.GRCh38.dna.primary_assembly.fa')
    #counts_fpath = join(cwd, 'counts', 'exon_skipping.filtered.total.csv')
    
    genome_fpath = sys.argv[1]
    counts_fpath = sys.argv[2]
    out_fpath = sys.argv[3]
    
    fastafile = Fastafile(genome_fpath)
    total = pd.read_csv(counts_fpath, index_col=0)
    features = get_5ss_sequences(total.index, fastafile)
    
    with open(out_fpath, 'w') as fhand:
        fieldnames = ['event_id', '5ss']
        writer = DictWriter(fhand, fieldnames=fieldnames)
        writer.writeheader()
        
        for record in tqdm(features, total=total.shape[0]):
            writer.writerow(record)

