"""
cellranger.py
=============

Helper functions for pipeline_cellranger.py

Code
====

"""

import pandas as pd
import numpy as np
import gzip
import shutil
import os
from cgatcore import pipeline as P


def get_stat(library_id, feature_type,
             samples,
             params):
    'return a statement to execute cellranger'

    fastq_path = samples.fastqs[[library_id, feature_type]]["path"]
    lib_dict = samples.libs["library_id"]

    feature_type = feature_type.to_lower()

    if feature_type == "vdj-t":
    
        cmd = '''cellranger vdj
                  --id vdj_t
                  --chain TR
                  --reference %(vdj_t_reference)
                  --fastqs %(fastq_path)s
               ''' % dict(PARAMS, **lib_dict, **locals())
    
    
        if params["vdj_t_r1-length"]:
            cmd += '''
                    --r1-length %(vdj_t_r1-length)s
                    '''
        if params["vdj_t_r2-length"]:
            cmd += '''
                    --r2-length %(vdj_t_r2-length)s
                    '''
        if params["vdj_t_inner-enrichment-primers"]:
            cmd += '''
                   --inner-enrichment-primers %(vdj_t_inner-enrichment-primers)s
                   '''
    
    elif feature_type == "vdj-b":
    
        cmd = '''cellranger vdj
                  --id vdj_b
                  --chain IG
                  --reference %(vdj_b_reference)
                  --fastqs %(fastq_path)s
               ''' % dict(PARAMS, **lib_dict, **locals())
    
    
        if params["vdj_b_r1-length"]:
            cmd += '''
                    --r1-length %(vdj_b_r1-length)s
                    '''
        if params["vdj_b_r2-length"]:
            cmd += '''
                    --r2-length %(vdj_b_r2-length)s
                    '''
        if params["vdj_b_inner-enrichment-primers"]:
            cmd += '''
                   --inner-enrichment-primers %(vdj_b_inner-enrichment-primers)s
                   '''
    
    else:
        raise ValueError("Unrecognised feature type")

    
    if params["cellranger_job_mode"] == "local":
    
        cmd += '''
               --jobmode local
               --localcores %(cellranger_localcores)s
               --localmem %(cellranger_localmem)s
        '''
    elif params["cellranger_job_mode"] == "cluster":
    
        cmd += '''
               --jobmode %(cellranger_job_template)s
               --maxjobs %(cellranger_max_jobs)s
               '''
        if params["cellranger_mempercore"]:
        
            cmd += '''
                   ---mempercore %(cellranger_mempercore)s
                   '''
    
    if feature_type.to_lower() in ['antibody_capture', 'gene_expression']:
    
        if params["cellranger_no-bam"]:
            cmd+= '''
                   --no-bam
                  '''
                  
        if params["cellranger_no-secondary"]:
            cmd+= '''
                   --no-secondary
                  '''
    
    cmd += '''
              --disable-ui
           '''

    return(cmd)

def get_counts(matrix_location, output_location,
               library_id):
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

    # define hash mapping of cellranger feature types
    # to three letter modality codes
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

            barcode_file = os.path.join(matrix_location,
                                        "barcodes.tsv.gz")

            barcode_outfile = os.path.join(out_path,
                                           "barcodes.tsv.gz")

            statement = '''zcat %(barcode_file)s
                           | awk 'BEGIN{FS="-"}{print $1"-%(library_id)s"}'
                           | gzip -c
                           >> %(barcode_outfile)s'''

            P.run(statement)

            out_mtx = os.path.join(out_path,
                                   "matrix.mtx")

            nrows = x.shape[0]
            ncol = input_mtx_dimensions[1]
            nentries = len(mtx_rows)

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

            # Gzip the outfile
            P.run('''gzip %(out_mtx)s''')

def contig_annotations(ctg_loc, out_loc, library_id):
    '''
    Fix the barcode syntax to "UMI-LIBRARY_ID" format
    '''

    outdir = os.path.dirname(out_loc)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    statement = '''cat %(ctg_loc)s
                   | sed '1s/barcode/barcode_id/g'
                   | sed 's/\(^[^-]*\)[^,]*\(.*$\)/\\1-%(library_id)s\\2/g'
                   | gzip -c
                   > %(out_loc)s
                '''

    P.run(statement)




def preprocess_cellranger_stats(infile, outfile):
    libraries = pd.read_csv(infile, sep='\t')
    libraries.set_index("library_id", inplace=True)

    ss = []

    for library_id in libraries.index:
        cellranger_summary = "/".join(["cellranger.multi.dir", library_id,
                                       "outs/per_sample_outs", library_id,
                                       "metrics_summary.csv"])
        map_sum = pd.read_csv(cellranger_summary, sep=',')

        # (need the group name to ensure entries are unique)
        map_sum = map_sum.replace(np.nan, 'na', regex=True)
        map_sum["id"] =  map_sum["Library or Sample"] + \
                         "_" + map_sum["Library Type"] + \
                         "_" + map_sum["Group Name"] + \
                         "_" + map_sum["Metric Name"]

        map_sum["id"] = [ x.replace(' ','_').lower() for x in map_sum["id"].values]
        map_sum.index = map_sum["id"]

        sub_map_sum = pd.DataFrame(map_sum["Metric Value"]).transpose()

        col_names = []
        for col in sub_map_sum.columns:

            x = sub_map_sum[col].values[0]

            if x[-1:]=="%":
                col_names.append(col + "_pct")
                x = x.replace("%","")
            else:
                col_names.append(col)

            x = x.replace(",","")


            sub_map_sum[col] = pd.to_numeric(x)
        sub_map_sum.columns = col_names
        sub_map_sum['library_id'] = library_id

        ss.append(sub_map_sum)

    summary = pd.concat(ss)
    summary.to_csv(outfile, sep = '\t', index=False)