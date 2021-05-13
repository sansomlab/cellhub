import pandas as pd
import numpy as np
import gzip
import shutil
import os
from cgatcore import pipeline as P

def get_counts(matrix_location, output_location):
    '''
    Post-process a cellranger multi count market matrix
    to generate individual per-modality matrices.

    This function assumes that different feature types
    are found in contiguous blocks of rows.
    '''

    if not os.path.exists(output_location):

        os.makedirs(output_location)

    feature_file = os.path.join(matrix_location, "features.tsv.gz")
    barcode_file = os.path.join(matrix_location, "barcode.tsv.gz")
    mtx_location = os.path.join(matrix_location, "matrix.mtx.gz")

    features = pd.read_csv(feature_file, sep="\t",
                           header=None, names=["id","name","type"])

    feature_types = set(features["type"])

    feature_names = {"Antibody Capture": "ADT",
                     "Gene Expression": "GEX",
                     "Multiplexing Capture": "HTO"}

    # TODO: if only one feature type skip the below and symlink instead.

    # get the row numbers from the market matrix
    # this is needed to determine which feature types are present

    fi = [] # fi = feature row index

    begin = True
    mtx_header = []
    with gzip.open(mtx_location, "rt", encoding="utf-8") as mtx:
        for line in mtx:
            if line.startswith("%"):
                mtx_header.append(line)
            else:
                if begin == True:
                    input_mtx_dimensions = [int(x) for x in line.split(" ")]
                    begin = False
                    continue
                fi.append(int(line.split(" ")[0]))


    fi = np.array(fi, dtype="int32")

    for feature_type in feature_types:

        if feature_type not in feature_names.keys():

            raise ValueError("Unrecognised feature type")

        x = features[features["type"]==feature_type]

        idx_start = min(x.index) + 1 # python is 0 indexed
        idx_end = max(x.index) + 1

        mtx_rows = fi[(fi >= idx_start) & (fi <= idx_end)]

        if(len(mtx_rows) > 0):

            # we have data for this feature type so we write it out

            out_path = os.path.join(output_location,
                                    feature_names[feature_type])

            if not os.path.exists(out_path):
                os.makedirs(out_path)

            # write out the features
            x.to_csv(os.path.join(out_path,
                                  "features.tsv.gz"),
                     sep="\t", header=None, index=False)

            # copy across the cell barcodes
            shutil.copyfile(os.path.join(matrix_location,
                                         "barcodes.tsv.gz"),
                            os.path.join(out_path,
                                         "barcodes.tsv.gz"))


            out_mtx = os.path.join(out_path,
                                   "matrix.mtx")

            nrows = x.shape[0]
            ncol = input_mtx_dimensions[1]
            nentries = len(fi)

            with open(out_mtx, "w") as mtx:

                mtx.write("".join(mtx_header))
                mtx.write(" ".join([str(nrows), str(ncol), str(nentries)]) + "\n")

            # subset the mtx matrix using AWK, because AWK is *fast*
            row_offset = idx_start - 1

            # append the matrix values, offsetting the column index
            statement = '''zcat %(mtx_location)s
                            | grep -v "^%%"
                            | awk 'NR>1 && $1>=%(idx_start)s && $1<=%(idx_end)s{
                                     print $1-%(row_offset)s,$2,$3}'
                            >> %(out_mtx)s;
                        '''
            P.run(statement)

            statement = '''gzip %(out_mtx)s'''
            P.run(statement)