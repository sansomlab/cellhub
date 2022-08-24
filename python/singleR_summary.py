import os
import re
import argparse
import anndata as ad
import pandas as pd
import logging
import sys
from pathlib import Path


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
parser.add_argument("--label_tables", default="none", type=str,
                    help="comma seperated list of label tables")
parser.add_argument("--outfile",default=1, type=str,
                    help="name of the output summary table")
args = parser.parse_args()

L.info("Running with arguments:")
print(args)

reference_abbreviations = {"HumanPrimaryCellAtlasData":"HPCA",
                "BlueprintEncodeData":"BlueEnc",
                "MouseRNAseqData":"MmRNAseq",
                "ImmGenData":"ImmGen",
                "DatabaseImmuneCellExpressionData": "ImmCell",
                "NovershternHematopoieticData":"NovHem",
                "MonacoImmuneData":"MonImm"}

tables = [x.strip() for x in args.label_tables.split(",")]

start = True
for tab in tables:
    
    ref = os.path.basename(Path(tab).parents[0])
    col_name = "singleR_" + reference_abbreviations[ref]
    
    labs = pd.read_csv(tab, sep="\t")
    labs.index=labs[["barcode","library_id"]]
    
    if start:
        out = labs[["barcode", "library_id", "pruned.labels"]]
        out.columns = ["barcode", "library_id", col_name]
        start = False
    else:
        out[col_name] = labs.loc[out.index, "pruned.labels"]

out.to_csv(args.outfile, sep="\t", index=False)
        

L.info("Complete")
