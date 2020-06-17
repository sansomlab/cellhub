import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import logging
import sys
import argparse
import pandas as pd

# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("scrublet")

# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--cellranger_dir", default="", type=str,
                    help="Matrix of counts")
parser.add_argument("--keep_barcodes_file", default=None, type=str,
                    help="")
parser.add_argument("--sample", default="sample", type=str,
                    help="sample or channel id")
parser.add_argument("--expected_doublet_rate", default=0.06, type=float,
                    help="the expected fraction of transcriptomes that are doublets, typically 0.05-0.1. Results are not particularly sensitive to this parameter.")
parser.add_argument("--min_counts", default=2, type=int,
                    help="")
parser.add_argument("--min_cells", default=3, type=int,
                    help="")
parser.add_argument("--min_gene_variability_pctl", default=85, type=int,
                    help="")
parser.add_argument("--n_prin_comps", default=30, type=int,
                    help="")
parser.add_argument("--outdir", default="", type=str,
                    help="output directory")
args = parser.parse_args()

# ########################################################################### #
# ############################# Run scrublet ################################ #
# ########################################################################### #

# Read in data

# Counts matrix
counts_matrix = args.cellranger_dir + "/matrix.mtx.gz"
counts_matrix = scipy.io.mmread(counts_matrix).T.tocsc()

# Features file
features_file = args.cellranger_dir + "/features.tsv.gz"
features = pd.read_table(features_file, compression='gzip', sep='\t', header=None)
features = features[1].to_numpy(dtype='<U19')

# Barcode file
barcode_file = args.cellranger_dir + "/barcodes.tsv.gz"
barcodes = pd.read_table(barcode_file, compression='gzip', sep='\t', header=None)

L.info('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
L.info('Number of features in gene list: {}'.format(len(features)))

# Subset barcodes and transform barcodes object to numpy array
keep_barcodes_file = args.keep_barcodes_file
if keep_barcodes_file != None:
  L.info('Keeping only barcodes in input file')
  keep_barcodes = pd.read_table(keep_barcodes_file, compression='gzip', sep='\t', header=None)
  L.info('Numer of barcodes to keep: {}'.format(keep_barcodes.shape[0]))
  idx = barcodes[0].isin(keep_barcodes[0])
  counts_matrix = counts_matrix[np.where(idx.values)[0],:]
  L.info('Final counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
  barcodes = barcodes.loc[np.where(idx.values)[0]][0].to_numpy()
  L.info('Final number of barcodes: {} '.format(len(barcodes)))
else:
  barcodes=barcodes[0].to_numpy()

# Initialize scrublet object
L.info('Initializing scrublet object with expected double rate {}'.format(args.expected_doublet_rate))
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=args.expected_doublet_rate)

# Run scrublet
L.info('Computing doublet scores and predicting doublets')
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=args.min_counts, 
                                                          min_cells=args.min_cells, 
                                                          min_gene_variability_pctl=args.min_gene_variability_pctl, 
                                                          n_prin_comps=args.n_prin_comps)

# Check whether a threshold for identifying doublets has been automatically identififed
if not hasattr(scrub, 'threshold_'):
  L.info('Failed to automatically identify doublet score threshold')
  L.info('Setting doublet score threshold to 1. WARNING: this is just a placeholder value')
  L.info('Values in scrub_predicted_doublets column in output will be set to NA')
  # Calling doublets
  scrub.call_doublets(threshold=1)
  # Creating output file with results  
  out = pd.DataFrame({'barcode': barcodes, 
                    'id': [args.sample] * len(barcodes),
                    'scrub_doublet_scores' : doublet_scores, 
                    'scrub_predicted_doublets' : "NA"})
else: 
  L.info('Calculated doublet score threshold: {}'.format(scrub.threshold_))
  # Creating output file with results  
  out = pd.DataFrame({'BARCODE': barcodes, 
                    'sample': [args.sample] * len(barcodes),
                    'scrub_doublet_scores' : doublet_scores, 
                    'scrub_predicted_doublets' : predicted_doublets})

L.info('Final output shape: {} rows, {} columns'.format(out.shape[0], out.shape[1]))

# Write output 
L.info('Writing file with doublet scores and predicted doublets')
out.to_csv(os.path.join(args.outdir, 
                        "_".join([args.sample, "scrublet.tsv.gz"])),
                        sep="\t", index = False)

# Save histogram
scrub.plot_histogram()
plt.savefig(os.path.join(args.outdir, 
                         "_".join([args.sample,"doublet_score_histogram.png"])))

# Run and plot UMAP
L.info('Running UMAP')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
L.info('Plotting UMAP')
scrub.plot_embedding('UMAP', order_points=True)
plt.savefig(os.path.join(args.outdir, 
                         "_".join([args.sample, "doublet_score_umap.png"])))

L.info("Complete")
