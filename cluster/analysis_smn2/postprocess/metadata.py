import os,sys
from typing import Optional, Dict, Tuple, List
import pandas as pd
import utils
import re

class metadata(object):
    """ The single class for reading the metadata Excel file.
    """
    __version__ = '0.0.1'
    def __init__(self, metadata_file):
        """ Initialize the metadata class
        
        Keyword Arguments:
        metadata_file {[str]} -- [description] (default: {'data/metadata.xlsx'})
        """
        self.metadata_file = metadata_file
    
    def file_existence(self):
        """ Check the metadata file is exist or not. If it is not exist, the code
        will be terminated.
        """
        file_exist = os.path.isfile(self.metadata_file)
        if file_exist:
            print('--> The file {} exist'.format(self.metadata_file))
        else:
            print('--> Error: Cant find metadata file {}'.format(self.metadata_file))
            sys.exit(1)

    def load_metadata(self):
        """ Load the Excel file in pandas framework. 
        This function will:
        - Load the excel file.
        - Get the sheet names in the excel file.
        - Drop rows in which the use=False in all the sheets.
        - Get the sample names
        """
        # Load the excel file
        self.excel_file = pd.read_excel(self.metadata_file, sheet_name=None)
        # Get the sheet names
        self.sheet_names = self.excel_file.keys()
        # Get LIDs
        self.LIDs = self.excel_file['LIDs']['LID']

    def validate_metadata(self):
        """ Step by step validation of metadata 
        before we start to proceed to next levels
        """
        print('--> Metadata Validation')

        # Check regx parse the sample sequences:
        reads_sheet = self.excel_file['reads']
        print('---> Check sample sequences and regex')
        for index, row in reads_sheet.iterrows():
            m = re.match(row['regex'], row['sequence'])
            if m:
                raw_list = m.groupdict()
                print('---> in read {} found the match {}'.format(row['read'], raw_list))
            if not m:
                print('---> in read {} regex is not correct'.format(row['read']))
                sys.exit(1)

