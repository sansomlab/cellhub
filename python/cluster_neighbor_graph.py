import os
import re
import argparse
import anndata as ad
import pandas as pd
import logging
import sys



# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

L = logging.getLogger(__name__)
log_handler = logging.StreamHandler(sys.stdout)
log_handler.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
log_handler.setLevel(logging.INFO)
L.addHandler(log_handler)
L.setLevel(logging.INFO)


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

L.info("parsing arguments")

parser = argparse.ArgumentParser()
parser.add_argument("--source_anndata", default="source.adata.h5ad", type=str,
                    help="path to the source anndata object")
parser.add_argument("--reduced_dims_name", default="pca", type=str,
                    help="the X_[name] of the reduced dimension")
parser.add_argument("--outfile",default=1, type=str,
                    help="name of the output anndata object")
parser.add_argument("--ncomps", default=20, type=int,
                    help="Number of reduced components to include in the knn computation")
parser.add_argument("--method", default="scanpy", type=str,
                    help="scanpy|hnsw|sklearn|precomputed (scanpy uses pynndescent)")
parser.add_argument("--k", default=20, type=int,
                    help="number of neighbors")
parser.add_argument("--metric", default="euclidean", type=str,
                    help="the distance metric")
parser.add_argument("--threads", default=4, type=int,
                    help="number of threads")
parser.add_argument("--fullspeed", default=False, action="store_true",
                    help="number of threads")

args = parser.parse_args()

L.info("Running with arguments:")
print(args)


# ########################################################################### #
# #########################  Make the anndata object ######################### #
# ########################################################################### #

sourceAdata = ad.read_h5ad(args.source_anndata, backed='r')

# make the anndata and populate the obs and vars
adata = ad.AnnData(shape = (sourceAdata.shape[0],0))
adata.obs.index = sourceAdata.obs.index.copy()

rdims = "X_" + args.reduced_dims_name

# The slot is renamed to ensure we avoid additional checks run by scvelo
# in the case that the rep name == "X_pca"
nn_rdims = rdims + "_copy_for_nn"
adata.obsm[nn_rdims] = sourceAdata.obsm[rdims][:,0:args.ncomps].copy()

# ########################################################################### #
# ####################### Nearest neighbor computation ###################### #
# ########################################################################### #

# Run neighbors
L.info( "Using " + str(args.k) + " neighbors")

if args.method == "scanpy":

    L.info("Computing neighbors using default scanpy method")
    from scanpy.preprocessing import neighbors

    neighbors(adata,
              n_neighbors = args.k,
              metric = args.metric,
              use_rep = nn_rdims)



elif args.method == "hnsw":

    L.info("Computing neighbors using hnswlib (with scvelo a la pegasus!)")
    # we use the neighbors function from scvelo (thanks!)
    # with parameters from pegasus (for a more exact result).

    from scvelo.pp import neighbors
    num_threads = (args.threads if args.fullspeed else 1)

    neighbors(adata,
              n_neighbors = args.k,
              n_pcs = None,
              use_rep = nn_rdims,
              knn = True,
              random_state = 0,
              method = 'hnsw',
              metric = args.metric,
              metric_kwds = {"M":20,
                             "ef":200,
                             "ef_construction":200},
              num_threads=num_threads)


elif args.method == "sklearn":

    L.info("Computing neighbors using sklearn")
    # we use the neighbors function from scvelo

    from scvelo.pp import neighbors

    neighbors(adata,
              n_neighbors = args.k,
              n_pcs = None,
              use_rep = nn_rdims,
              knn = True,
              random_state = 0,
              method = 'sklearn',
              metric = args.metric,
              # metric_kwds = {}
              num_threads=args.threads)


elif args.method == "precomputed":

   L.info("Using precomputed neighbors")

   #add knn to adata
   
   adata.uns["neighbors"] = sourceAdata.uns["neighbors"]
   adata.obsp = sourceAdata.obsp
   
   #rename rep_use to 'X___'
   from scanpy._utils import NeighborsView
   neighbors = NeighborsView(adata, "neighbors")
        
   neighbors["params"]["use_rep"]='X_'+args.reduced_dims_name+'_copy_for_nn'
    
   adata.uns["neighbors"]["params"]=neighbors["params"] 

else:
    raise ValueError("nn method not recognised")

# compute clusters
adata.write(args.outfile)                         

L.info("Complete")
