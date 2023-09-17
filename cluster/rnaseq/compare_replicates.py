import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from itertools import combinations
from scipy.stats import pearsonr

fpath = 'counts/exon_skipping.filtered.psi.csv'
psi = pd.read_csv(fpath, index_col=0)

fig, subplots = plt.subplots(5, 6, figsize=(12, 10), sharex=True, sharey=True)
subplots = subplots.flatten()

reps = list(range(1, 6))
conditions = ['dmso', 'ris', 'bran']

rhos = []
k = 0
for condition in conditions:
    for i, j in combinations(reps, 2):
        sample1, sample2 = '{}_{}'.format(condition, i), '{}_{}'.format(condition, j)
        p = psi[[sample1, sample2]].dropna().values
        x, y = p[:, 0], p[:, 1]
        r1 = pearsonr(x, y)[0]
        
        p = p[np.all((p > 0.05) & (p < 0.95), axis=1),:] * 100
        x, y = p[:, 0], p[:, 1]
        r2 = pearsonr(x, y)[0]
        rhos.append({'sample1': sample1, 'sample2': sample2, 'r': r1, 'r_alternative': r2})
        
        axes = subplots[k]
        sns.histplot(x=x, y=y, cmap='viridis', ax=axes)
        axes.set(xlim=(0, 100), ylim=(0, 100))
        k += 1
        
rhos = pd.DataFrame(rhos)

fig.supxlabel('PSI Sample 1', ha='center', va='center')
fig.supylabel('PSI Sample 2', ha='center', va='center')
fig.tight_layout()
fig.savefig('counts/samples_psi.png')
rhos.to_csv('counts/samples_psi.csv', index=False)

print('All exons average Pearson r={:.3f}; SD={:.3f}'.format(rhos['r'].mean(), rhos['r'].std()))
print('Exons with 5<PSI<95 average Pearson for  r={:.3f}; SD={:.3f}'.format(rhos['r_alternative'].mean(), rhos['r_alternative'].std()))
