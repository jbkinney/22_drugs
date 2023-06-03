import os 
import pandas as pd
from metadata import *
import numpy as np
from datetime import date
today = date.today()
dir_path = os.path.dirname(os.path.realpath(__file__))

metadata_file    = dir_path + '/../metadata_smn2.xlsx'
count_file_dir   = dir_path + '/../intermediate/counts'
report_file_dir  = dir_path + '/../intermediate/reports'
cipher_dir       = dir_path + '/cipher'
results_dir      = dir_path + '/results'
psi_dir          = dir_path + '/../psi_data'
report_file      = dir_path + '/report.csv'
