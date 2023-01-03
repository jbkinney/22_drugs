import sys
import os
import pandas as pd
import numpy as np
import re
import glob
from pandas.api.types import is_string_dtype
from typing import Optional, Dict, Tuple, List
import time

# User options
illumina_run_dir    = '/sonas-hs/kinney/hpc/home/jkinney/big_data/illumina_runs'
metadata_file       = './metadata_smn2_cipher.xlsx'
tmp_dir             = 'tmp'
interm_dir          = 'intermediate'
reads_per_split     = int(1E6)
clean_intermediate  = False
make_split_files    = True
make_features_files = True
make_counts_files   = True
merge_feature       = False
clean_tmp = True
use_cluster = True
