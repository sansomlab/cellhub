import scanpy as sc
import pandas as pd
import os
import yaml
import argparse
import h5py
import anndata

parser = argparse.ArgumentParser()
parser.add_argument("--task-yml", default="", type=str,
                    help="yml for this task")
args = parser.parse_args()


def addMetadata(adata, metadata_infile, id_col):
    metadata_add = pd.read_csv(metadata_infile, sep = "\t",
                               compression="gzip")
    metadata_add.drop_duplicates(inplace=True)
    new_obs = pd.merge(adata.obs.copy(), metadata_add,
                       on= id_col, how='left')
    new_obs.set_index('barcode_id', inplace=True)
    new_obs = new_obs.loc[adata.obs.index,:]
    adata.obs = new_obs
    return adata


##### read inputs #####

# Read YAML file
with open(args.task_yml, 'r') as stream:
    opt = yaml.safe_load(stream)

print("Running with arguments:")
print(opt)

if not os.path.exists(opt["outdir"]):
    os.makedirs(opt["outdir"])

print(opt["infile_h5"])
adata = sc.read_h5ad(opt["infile_h5"])

adata.obs['barcode_id'] = adata.obs.index.values
if 'metadata_file' in opt.keys():
    print("adding metadata")
    adata = addMetadata(adata = adata,
                        metadata_infile = opt['metadata_file'],
                        id_col = opt['metadata_id'])
print(adata.obs.columns)

#### get necessary data from regressed+integrated object
scaled_data = pd.DataFrame(adata.X[:,adata.var["highly_variable"]].copy())
barcodes = adata.obs.index
genes = adata.var.index[adata.var["highly_variable"]].copy()

## save harmony components
dim = "X_" + opt["dim_name"]
comp_out = adata.obsm[dim].copy()
#harmony_out.columns = ["harmony_" + str(i) for i in range(1,harmony_out.shape[1]+1)]

metadata_out = adata.obs.copy()
# if barcode exists, this is likely incorrect -> use index instead
if 'barcode' in metadata_out.columns:
    del metadata_out['barcode']
metadata_out.index.name = 'barcode'
metadata_out.reset_index(inplace=True)
metadata_out.to_csv(os.path.join(opt["outdir"], "metadata.tsv.gz"),
                    sep="\t", index=False, compression="gzip")

print("create h5 file")
hf = h5py.File(os.path.join(opt["outdir"], "scaled.h5"), 'w')
if not "seurat_data_only" in opt.keys() :
    print("add scaled data to h5file")
    hf.create_dataset('scaled_data', data=scaled_data)
hf.create_dataset('barcodes', data=barcodes)
hf.create_dataset('genes', data=genes)
hf.create_dataset('comp', data=comp_out)

hf.close()

### recreate the full object to match h5ad requirements
## counts need to be in raw
if 'infile_full' in opt.keys():
    print("Create new anndata object")
    outfile_anndata = os.path.join(opt["outdir"], "full_anndata.h5ad")
    adata_full = sc.read_h5ad(opt["infile_full"])
    adata_new = anndata.AnnData(X=adata_full.layers['counts'].copy(),
                                obs=adata_full.obs.copy(),
                                var=adata_full.var.copy())
    adata_new.raw = adata_new
    adata_new.obs['barcode_id'] = adata_new.obs.index.values
    if 'metadata_file' in opt.keys():
        print("adding metadata")
        adata_new = addMetadata(adata = adata_new,
                            metadata_infile = opt['metadata_file'],
                            id_col = opt['metadata_id'])
    print("New anndata object has the following obs: ")
    print(adata_new.obs.columns)
    adata_new.X = adata_full.X
    adata_new.write(outfile_anndata)

print("Completed")
