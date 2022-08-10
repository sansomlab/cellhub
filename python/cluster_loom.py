import os
import re
import argparse
import anndata as ad
import loompy
import pandas as pd
import numpy as np
import scipy as sc
import logging
import sys

# <------------------------------ Logging ------------------------------------>

L = logging.getLogger(__name__)
log_handler = logging.StreamHandler(sys.stdout)
log_handler.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
log_handler.setLevel(logging.INFO)
L.addHandler(log_handler)
L.setLevel(logging.INFO)

# <------------------------------ Arguments ---------------------------------->

L.info("parsing arguments")

parser = argparse.ArgumentParser()
parser.add_argument("--anndata", default="adata.h5ad", type=str,
                    help="path to the anndata object")
parser.add_argument("--layer", default="adata.h5ad", type=str,
                    help='The layer to export, either "X" or the name of a layer')
parser.add_argument("--loomdir", default="X_pca", type=str,
                    help="path to save the loom object")

args = parser.parse_args()

L.info("Running with arguments:")
print(args)

# <------------------------------ Functions ---------------------------------->

def _chunked_loom_export(source_anndata = None, 
                         loom_file = None, 
                         anndata_layer = "X", 
                         chunk_size=10000):
    '''
    Perform a chunked export of a single matrix from an anndata file to a loom file.
    '''
    

    loom_conn = loompy.new(loom_file)

    ncells = source_anndata.shape[0]
    nchunks = int(np.ceil(ncells/chunk_size))
        
    row_attrs = {"gene_name": np.array(source_anndata.var.index)} 
    
    for chunk in range(0, nchunks):
    
        L.info("Processing chunk: " + str(chunk))
        
        chunk_start = chunk * chunk_size
        
        if chunk == (nchunks - 1):
            chunk_end = ncells
        else:
            chunk_end = ((chunk + 1) * chunk_size )
            
        L.info("... chunk start: " + str(chunk_start) + ", chunk end:" + str(chunk_end))
        
        if anndata_layer == "X":
            data = source_anndata.X[chunk_start:chunk_end,:]
        else:
            data = source_anndata.layers[anndata_layer][chunk_start:chunk_end,:]
                        
        if isinstance(data, np.ndarray):
            data = np.transpose(data)
            
        elif isinstance(data, sc.sparse.csr.csr_matrix):
            data = data.transpose()
            data = data.toarray()
            
        else:
            L.info("matrix type: ")
            print(type(data))
            
            raise ValueError("data type not recognised")
        
        col_attrs = {"barcode_id": np.array(source_anndata.obs.barcode_id[chunk_start:chunk_end])}
        
        print("adding data slice to loom")
        if chunk == 0:
            loom_conn.add_columns(data, col_attrs=col_attrs, row_attrs=row_attrs)
          
        else:
            loom_conn.add_columns(data, col_attrs=col_attrs)
            
    loom_conn.close()


# <--------------------------- Make the loom(s) ------------------------------>

L.info("Reading the anndata: " + args.anndata)
source_adata = ad.read_h5ad(args.anndata, backed='r') 

if args.layer == "X":
    if source_adata.X is None:
        raise ValueError("adata.X does not exist")
    
else:
    if source_adata.layers[args.layer] is None:
        raise ValueError("The specificed layer does not exist")


loom_file = os.path.join(args.loomdir, args.layer + ".loom")  
L.info('Exporting the "' + args.layer + '" layer to: ' + loom_file)

_chunked_loom_export(source_adata, 
                    loom_file = loom_file, 
                    anndata_layer=args.layer, 
                    chunk_size=20000)

L.info("complete")