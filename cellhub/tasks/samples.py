'''
samples.py
======

Overview
--------


- read in samples to pandas table
- read in fastqs to pandas table

- make sample object
     - has metadata
     - has seqdata
        -  feature type
              - path [can be appended]

Usage
-----

Registering outputs on the service endpoint
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




Class and method documentation
------------------------------

'''

import yaml
import os
import shutil
import re
import copy
import pandas as pd
from pprint import pprint

def check_cols(pd_frame, req_columns_list,
               table_name="default"):

    for col_name in req_columns_list:
        if col_name not in pd_frame.columns:
            raise ValueError("required column: '" + col_name "' missing "
                             "in " + table_name + " table")


class samples():
    '''
    A class for defining samples.

    When initialising an instance of the class, the pipeline name
    is passed e.g.::

      x = cellhub.tasks.api.register("cell_qc")

    .. note:: pipeline names are sanitised to replace spaces, underscores and hypens with periods.
    '''

    def __init__(self, pipeline = None, fastq_tsv = None, library_tsv = None):

        if pipeline is None or pipeline == "":
            raise ValueError("a pipeline name must be specified")

        self.pipeline = re.sub("[ \-_]",".",pipeline)

        libs = pd.read_csv(library_tsv, sep="\t")
        lib_cols = ["sample_name", "chemistry"]
        
        check_cols(libs,lib_cols,
                   "libraries.tsv")
        
        fastqs = pd.read_csv(fastq_tsv)
        
        check_cols(fastqs, ["path","library_id","feature_type","lane"],
                   "fastqs.tsv")
        
        if len(libs["library_id"].values) > len(set(meta["library_id"].values)):
            raise ValueError("Non-unique library_ids provided")
            
        libraries = libs["library_id"].values
        
        libs.index = libs["library_id"]
        self.libs = libs.to_dict(orient="index")
        
        fastqs.sort_values(["library_id",
                    "feature_type",
                    "lane"], inplace=True)

        f = (fastqs.groupby(["libarary_id","feature_type"])
                .agg({'library_id':'first', 'feature_type':'first', 
                      'path':','.join}))
        
        f.index = f["library_id"]
    
        self.fastqs = f


        



