import pandas as pd
from metadata import *
import numpy as np
from datetime import date
today = date.today()


metadata_file    = '/grid/kinney/home/kooshkb/pipeline_mpsa/ikbkap_drug_110321/ikbkap_drug_110321.xlsx'
count_file_dir   = '/grid/kinney/home/kooshkb/pipeline_mpsa/ikbkap_drug_110321/intermediate/counts/'
report_file_dir  = '/grid/kinney/home/kooshkb/pipeline_mpsa/ikbkap_drug_110321/intermediate/reports'
cipher_dir     = './cipher'
results_dir    = './results_'+str(today)
report_file    = '../report_'+str(today)+'.csv'
