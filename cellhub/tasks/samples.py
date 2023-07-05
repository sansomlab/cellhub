'''
samples.py
======

Overview
--------


- read in samples to pandas table
- read in libraries to pandas table

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


class samples():
    '''
    A class for defining the samples and libraries present in
    a 10x experiment.
    '''

    def __init__(self,
                 sample_tsv = None, 
                 library_tsv = None):

        samples = pd.read_csv(sample_tsv, sep="\t")
        sample_cols = ["sample_id", "library_id", 
                       "chemistry", "expect_cells"]
        
        check_cols(samples,sample_cols,
                   "samples.tsv")
        
        libs = pd.read_csv(library_tsv, sep="\t")
        
        check_cols(libs, ["library_id","feature_type",
                            "sample","fastq_path"],
                   "libraries.tsv")
        
        self.feature_types = ["Gene Expression",
                         "Antibody Capture",
                         "CRISPR Guide Capture",
                         "Custom"]
        
        self.vdj_types = ["VDJ-T", "VDJ-B"]
        
        check_values(libs, "feature_type", 
                     self.feature_types + self.vdj_types)
        
        if len(samples["library_id"].values) > len(set(libs["library_id"].values)):
            raise ValueError("Non-unique library_ids provided")
            
        libraries = samples["library_id"].values
        
        samples.index = samples["library_id"]
        self.samples = samples.to_dict(orient="index")
        
        libs.sort_values(["library_id",
                    "feature_type",
                    "sample"], inplace=True)
        
        self.libs = libs

    def get_samples_and_fastqs(self, library_id, feature_type):
        '''
        Return sample(s) and path(s) to fastq(s) for given library and
        feature type. Used to build job statements for cellranger vdj.
        
        If multiple sequence files are provide Samples and fastqc_paths
        will be merged into comma seperated lists, see:
        https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/vdj
        
        '''
        
        agg_libs = (self.libs.groupby(["library_id","feature_type"])
                .agg({'library_id':'first', 'feature_type':'first', 
                      'sample':','.join, 'fastq_path':','.join}))
        
        agg_libs.index = agg_libs["library_id"]
                
        x = agg_libs[(agg_libs["library_id"]==library_id) &
                        (agg_libs["feature_type"]==feature_type)]
        
        if x.shape[0] > 1:
            print(x)
            print(x.shape)
            raise ValueError("Sample and path aggregation failed, check input files.")
        
        fastq_path = x["fastq_path"].values[0]
        # check that the fastq paths are different (see note in pipeline_cellranger.py)
        fqps = [x.strip() for x in fastq_path.split(",")]
        
        if len(set(fqps)) < len(fqps):
        
            raise ValueError("Duplicate FASTQ paths detected for data from different"
                             " flow cells. This is not supported by the 'cellranger vdj' command." 
                             " VDJ data from different flow cells must be arranged in different" 
                             " folders. See note in pipeline_cellranger.py")

        sample = x["sample"].values[0]
        
        return({"fastq_path":fastq_path, "sample":sample})
        

    def feature_barcode_libraries(self):
        '''
        Return a list of libraries with Gene Expression/Antibody
        Capture data
        '''

        return(set(
            self.libs[self.libs["feature_type"].isin(
            self.feature_types)]["library_id"].values
        ))
        
    def vdj_libraries(self):
        '''
        Return a list of libraries with VDJ data
        '''
        return(set(
            self.libs[self.libs["feature_type"].isin(
            self.vdj_types)]["library_id"].values
        ))
        
    def vdj_t_libraries(self):
        '''
        Return a list of libraries with VDJ data
        '''
        return(set(
            self.libs[self.libs["feature_type"]=="VDJ-T"]["library_id"].values
        ))
        
    def vdj_b_libraries(self):
        '''
        Return a list of libraries with VDJ data
        '''
        return(set(
            self.libs[self.libs["feature_type"]=="VDJ-B"]["library_id"].values
        ))
        
    def lib_types(self, library_id):
    
        return(set([x for x in 
                    self.libs[self.libs["library_id"]==library_id]["feature_type"].values
                    ]))
        

    def write_csv(self, library_id, outfile_path):
        '''Return a Libraries CSV file for pipeline cellranger count
           in format:
           
           fastqs,sample,library_type
            /opt/foo/,GEX_sample1,Gene Expression
            /opt/bar/,GEX_sample1,Gene Expression
            /opt/foo/,Ab_sample1,Antibody Capture
            /opt/foo/,CRISPR_sample1,CRISPR Guide Capture
            
            see: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis#libraries-csv
        '''
        
        out = self.libs[self.libs["library_id"]==library_id]
        
        out = out[out["feature_type"].isin(self.feature_types)]
        
        out = out[["fastq_path", "sample","feature_type"]]
        out.columns = ["fastqs","sample","library_type"]
        
        out.to_csv(outfile_path, index=False, sep=",")
        
    #def library_parameters(self, library_id):
    
        
    
            

    
    
        



