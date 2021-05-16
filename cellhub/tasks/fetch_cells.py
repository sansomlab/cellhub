import os
import sys
import shutil
import gzip
import hashlib
from pathlib import Path
from cgatcore import pipeline as P

# adapted from:
# https://stackoverflow.com/questions/3431825/
# generating-an-md5-checksum-of-a-file

def md5gz(fname):

    hash_md5 = hashlib.md5()
    with gzip.open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def get_cell_subset(barcodes,
                    modality,
                    matrix_id,
                    outdir,
                    PARAMS,
                    data_subset="filtered",):
    '''
    subset a mtx matrix
    '''
    #api_location=os.path.join(PARAMS["cellhub_location"],
    #                          "api"),

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if ("G" in PARAMS["resources_memory"] or
        "M" in PARAMS["resources_memory"] ):
        job_memory = PARAMS["resources_memory"]

    log_file = os.path.join(outdir, "extract_cells.log")

    to_cluster = True

    source_matrix_dir = os.path.join(PARAMS["cellhub_location"],
                                     "api",
                                     "cellranger.multi",
                                     modality,
                                     data_subset,
                                     matrix_id,
                                     "mtx")

    statement = '''Rscript %(cellhub_code_dir)s/R/extract_cells.R
                   --cells=%(barcodes)s
                   --matrixdir=%(source_matrix_dir)s
                   --matrixid=%(matrix_id)s
                   --outdir=%(outdir)s
                   &> %(log_file)s
                '''

    P.run(statement)


def merge_subsets(infiles, outfile, PARAMS):
    '''merge a set of market matrix files'''

    ncells = 0

    # get the dimensions of all the market matrix files
    mtx_specs = {}

    mtx_encoding = "us-ascii"

    for subset in infiles:

        matrix_id = os.path.basename(Path(subset).parent)

        with gzip.open(subset, "r") as mtx:
            for i, line in enumerate(mtx, 1):
                if i == 1:
                    mtx_header = line  # .decode("us-ascii").strip()
                if i == 2:
                    line_str = line.decode(mtx_encoding)
                    nrow, ncol, nnonzero = line_str.strip().split(" ")
                    mtx_specs[matrix_id] = {"nrow": int(nrow),
                                            "ncol": int(ncol),
                                            "nnonzero": int(nnonzero),
                                            "mtx_file": subset}
                    break

    matrices_to_merge = tuple(mtx_specs.keys())

    mtx_outfile = outfile[:-len(".gz")]
    barcodes_outfile = os.path.join(os.path.dirname(outfile),
                                    "barcodes.tsv")

    features_outfile = os.path.join(os.path.dirname(outfile),
                                    "features.tsv.gz")

    if os.path.exists(mtx_outfile) or \
       os.path.exists(mtx_outfile + ".gz"):
        raise ValueError("mtx outfile already exists")

    if os.path.exists(barcodes_outfile) or \
       os.path.exists(barcodes_outfile + ".gz"):
        raise ValueError("barcodes outfile already exists")

    if os.path.exists(features_outfile):
        raise ValueError("features outfile already exists")

    # construct the header of the market matrix file.
    nrows = []
    total_ncol = 0
    total_nnonzero = 0

    for matrix_id in matrices_to_merge:
        nrows.append(mtx_specs[matrix_id]["nrow"])
        total_ncol += mtx_specs[matrix_id]["ncol"]
        total_nnonzero += mtx_specs[matrix_id]["nnonzero"]

    if len(set(nrows)) > 1:
        raise ValueError("the input matrices have different numbers of rows!")

    out_spec = " ".join([str(nrows[0]),
                         str(total_ncol),
                         str(total_nnonzero)]) + "\n"

    # write the header of the market matrix file
    with open(mtx_outfile, "wb") as out:
        out.write(mtx_header)
        out.write(out_spec.encode(mtx_encoding))

    column_offset = 0

    feature_file_checksums = []

    statement = ""

    for matrix_id in matrices_to_merge:

        mtx_file = mtx_specs[matrix_id]["mtx_file"]

        barcodes_file = os.path.join(os.path.dirname(mtx_file),
                                     "barcodes.tsv.gz")

        # append the matrix values, offsetting the column index
        statement += '''zcat %(mtx_file)s
                       | awk 'NR>2{print $1,$2 + %(column_offset)s,$3}'
                       >> %(mtx_outfile)s;
                    ''' % locals()

        # append the barcodes, adding the matrix identifier
        # | awk '{print $1"-1-%(matrix_id)s"}'
        # SNS may 2021: no longer necessary as we switch to
        # using the barcode_id which is set upstream in
        # cellranger_multi.
        # (TODO additional work may be needed to guarentee uniqueness
        # between different cellhub instances).
        statement += '''zcat %(barcodes_file)s
                       >> %(barcodes_outfile)s;
                    ''' % locals()

        # increase the column offset by the number of appended columns
        column_offset += mtx_specs[matrix_id]["ncol"]

        features_file = os.path.join(os.path.dirname(mtx_file),
                                     "features.tsv.gz")

        feature_file_checksums.append(md5gz(features_file))

    if ("G" in PARAMS["resources_memory"] or
        "M" in PARAMS["resources_memory"] ):
        job_memory = PARAMS["resources_memory"]

    # run the job.
    P.run(statement)

    # check that all of the subsets have identical features.
    if len(set(feature_file_checksums)) != 1:
        raise ValueError("The matrices have different features")
    else:
        print("The matrices have the same features")

    if ("G" in PARAMS["resources_memory"] or
        "M" in PARAMS["resources_memory"] ):
        job_memory = PARAMS["resources_memory"]

    # compress the outfiles
    statement = '''gzip %(mtx_outfile)s;
                   gzip %(barcodes_outfile)s
                '''

    P.run(statement)

    # copy over the features to the output directory
    shutil.copyfile(features_file, features_outfile)


def downsample_gex(agg_matrix_dir, cell_metadata,
                   outdir, PARAMS):

    group_var = PARAMS["GEX_downsample_group_var"]
    downsample_function= PARAMS["GEX_downsample_function"]

    log_file = os.path.join(outdir,
                            "downSampleTranscripts.log")

    if ("G" in PARAMS["resources_memory"] or
        "M" in PARAMS["resources_memory"] ):
        job_memory = PARAMS["resources_memory"]

    statement = '''Rscript %(cellhub_code_dir)s/R/downSampleTranscripts.R
                   --genexpdir=%(agg_matrix_dir)s
                   --cellmetadata=%(cell_metadata)s
                   --cellgroupvar=%(group_var)s
                   --downsample=%(downsample_function)s
                   --outdir=%(outdir)s
                   &> %(log_file)s
                '''

    P.run(statement)



def export_anndata(mtx_dir,
                   obs_file,
                   outdir,
                   matrix_name,
                   PARAMS):

    if not os.path.exists(outdir):

        os.mkdir(outdir)

    log_file = os.path.join(outdir,
                            matrix_name  + ".log")

    if ("G" in PARAMS["resources_memory"] or
        "M" in PARAMS["resources_memory"] ):
        job_memory = PARAMS["resources_memory"]

    job_threads = 1

    statement = '''python %(cellhub_code_dir)s/python/convert_mm_to_h5ad.py
                          --mtxdir10x=%(mtx_dir)s
                          --obsdata=%(obs_file)s
                          --obstotals=total_UMI
                          --outdir=%(outdir)s
                          --matrixname=%(matrix_name)s
                    > %(log_file)s
                 '''

    P.run(statement)
