'''
dehash.py
=========

Helper functions for pipeline_dehash.py

code
====

'''

import pandas as pd
import numpy as np
import os

def parse_gmmdemux(gmm_dir,
                   results_dir,
                   library_id,
                   hto_names):
    '''
    Parse the output of GMM-demux
    '''

    # index = barcode_id
    x = pd.read_csv(os.path.join(gmm_dir,
                                 "simple",
                                 "GMM_simplified.csv"),
                    index_col=0)

    x.columns = [ "gmm_" + xc.lower() for xc in x.columns]

    # index = Cluster_id
    y = pd.read_csv(os.path.join(gmm_dir,
                                 "simple",
                                 "GMM_simplified.config"),
                    index_col=0, names=["call"])

    y.columns = [ "gmm_" + yc.lower() for yc in y.columns]

    # clean up the white space ...
    y["gmm_call"] = y["gmm_call"].str.strip()

    out = pd.merge(x, y,
                   how="left",
                   left_on="gmm_cluster_id",
                   right_index=True)

    out["barcode"] = out.index
    out["library_id"] = library_id

    out["gmm_singlet"] = out["gmm_call"].isin(hto_names)

    out.to_csv(os.path.join(results_dir,
               library_id +".tsv.gz"),
               sep="\t",
               index=False)
