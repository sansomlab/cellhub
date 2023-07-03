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
            raise ValueError("required column: '" + col_name + "' missing "
                             "in " + table_name + " table")


def check_values(pd_frame, col, allowed):
    '''utility function for sanity checking columns'''
    if not all([x in allowed for x in pd_frame[col].values]):
        raise ValueError("Only the following values are allowed in column '"
                         + col + "': " + ",".join(allowed))


class lib():
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
        lib_cols = ["library_id", "chemistry", "expect_cells"]
        
        check_cols(libs,lib_cols,
                   "libraries.tsv")
        
        fastqs = pd.read_csv(fastq_tsv, sep="\t")
        print(fastqs)
        
        check_cols(fastqs, ["library_id","feature_type",
                            "seq_lane","fastq_path"],
                   "fastqs.tsv")
        
        self.feature_types = ["Gene Expression",
                         "Antibody Capture",
                         "CRISPR Guide Capture",
                         "Custom"]
        
        self.vdj_types = ["VDJ-T", "VDJ-B"]
        
        check_values(fastqs, "feature_type", 
                     self.feature_types + self.vdj_types)
        
        if len(libs["library_id"].values) > len(set(fastqs["library_id"].values)):
            raise ValueError("Non-unique library_ids provided")
            
        libraries = libs["library_id"].values
        
        libs.index = libs["library_id"]
        self.libs = libs.to_dict(orient="index")
        
        fastqs.sort_values(["library_id",
                    "feature_type",
                    "seq_lane"], inplace=True)

        f = (fastqs.groupby(["libarary_id","feature_type"])
                .agg({'library_id':'first', 'feature_type':'first', 
                      'fastq_path':','.join}))
        
        f.index = f["library_id"]
        
        self.fastqs = f

    def feature_barcode_libraries(self):
        '''
        Return a list of libraries with Gene Expression/Antibody
        Capture data
        '''

        return(set(
            self.fastqs[self.fastqs["feature_type"].isin(
            self.feature_types)]["library_id"].values
        ))
        
    def lib_types(self, library_id):
    
        return(set([x for x in 
                    self.fastqs[self.fastqs["library_id"]==library_id]["feature_type"].values
                    ]))
        

    def write_csv(self, library_id, outfile_path):
        '''Return a Libraries CSV file for pipeline cellranger count
           in format:
           
           fastqs,sample,library_type
            /opt/foo/,GEX_sample1,Gene Expression
            /opt/bar/,GEX_sample1,Gene Expression
            /opt/foo/,Ab_sample1,Antibody Capture
            /opt/foo/,CRISPR_sample1,CRISPR Guide Capture
        '''
        
        out = self.fastqs[self.fastqs["library_id"]==library_id]
        out = out[["fastq_path", "library_id","feature_type"]]
        out.columns = ["fastqs","sample","library_type"]
        
        out.to_csv(outfile_path, index=False, sep=",")
        
    #def library_parameters(self, library_id):
    
        
    
            

    
    
        



