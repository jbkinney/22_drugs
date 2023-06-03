#!/usr/bin/env python
from os.path import abspath, join, dirname

import sys
import pandas as pd

basedir = abspath(dirname(__file__))

prefix = join(basedir, 'bam')
suffix = 'Aligned.sortedByCoord.out.bam'
fpath = join(basedir, 'design.csv')

design = pd.read_csv(fpath, index_col=0)

group = sys.argv[1]

if group == 'control':
    samples = design.index[design[['Risdiplam', 'Branaplam']].sum(1) == 0].values
elif group == 'treated':
    samples = design.index[design[['Risdiplam', 'Branaplam']].sum(1) > 0].values
else:
    raise ValueError('Group not allowed {}'.format(group))

line = ','.join([join(prefix, '{}{}'.format(sample, suffix)) for sample in samples])

fpath = join(basedir, '{}.txt'.format(group))
with open(fpath, 'w') as fhand:
    fhand.write(line)

