#!/usr/bin/env python

import torch
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import cell2location
import pandas as pd
import os
import subprocess
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
import argparse
parser = argparse.ArgumentParser()


def int_or_none(value):
    if value.lower() == 'none':
        return None
    try:
        return int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid int value: '{value}'")


def print_gpu_usage():
    result = subprocess.run(['nvidia-smi', '--query-gpu=memory.used,memory.free', '--format=csv,nounits,noheader'], 
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE, 
                            check=True, 
                            text=True)
    print(result.stdout)

####figure out which ones we need, since we are setting them up in Main ####
parser.add_argument("--project_folder", type=str, help="Indicate the path to the project")
parser.add_argument("--out_dir", type=str, help="output_subdirectory in processed")
parser.add_argument("--sc_path", type=str, help="Indicate the path with the scRNASeq object")
parser.add_argument("--sc_batch", type=str, help="Indicate the column in the metadata of the scRNASeq object with the sample id to remove batch effect")
parser.add_argument("--sc_covariates", type=str, help="other variables that should be included to remove technical effects")
parser.add_argument("--sc_annotation_column", type=str, help="Indicate the column in the metadata of the scRNASeq object with the annotation")
parser.add_argument("--st_path", type=str, help="Indicate the path with the spatial h5ad object")
parser.add_argument("--st_batch", type=str, help="batch key in spatial data")
parser.add_argument("--st_covariates", type=str, help="other variables that should be included to remove batch effects")
parser.add_argument("--sc_batch_size", type=int, default=2048,  help="batch size for reference model training")
parser.add_argument("--st_batch_size", type=int_or_none, default=None,  help="batch size for spatial model training; default is to train on all spots")
parser.add_argument("--epochs_ref", type=int, default=500,  help="epochs for training the reference single cell model")
parser.add_argument("--epochs_st", type=int, default=20000,  help="epochs for training the spatial model")

args = parser.parse_args()
print_gpu_usage()

# SET VARIABLES
output_folder = f"{args.project_folder}/processed/{args.out_dir}/"
reports_folder = f'{args.project_folder}/reports/{args.out_dir}' #plots

# create output directories if they don't exist
os.makedirs(f"{output_folder}", exist_ok=True)
os.makedirs(f"{reports_folder}", exist_ok=True)

# where to save models
ref_model = f'{output_folder}/reference_signatures'
st_model = f'{output_folder}/cell2location_map'

# read sc DATA & input
print("Reading in single cell data")

sc_path = args.sc_path
sc_annotation_column = args.sc_annotation_column
sc_batch = args.sc_batch

covariate_list = args.sc_covariates.split(",")
sc_categorical_covariate = [x.strip() for x in covariate_list]
if len(sc_categorical_covariate) == 0:
    sc_categorical_covariate = None

st_batch = args.st_batch
st_covariate_list = args.st_covariates.split(",")
st_categorical_covariate = [x.strip() for x in st_covariate_list]
if len(st_categorical_covariate) == 0:
    st_categorical_covariate = None

adata_sc = sc.read_h5ad(sc_path)

# read ST 
print("Reading in spatial data")

ST_adata = sc.read_h5ad(args.st_path)
print("GPU usage after loading ST and SC data")
print_gpu_usage()

# finish setting up single cell and spatial anndata objects
adata_sc.var['SYMBOL'] = adata_sc.var_names # feature ids are gene symbols

# make sure raw counts are in single cell X slot 
if "counts" in adata_sc.layers:
    adata_sc.X = adata_sc.layers["counts"]
else:
    # assume adata_sc.X is the raw counts
    pass

# ST
ST_adata.var['SYMBOL'] = ST_adata.var_names #set gene as symbol
# not using ensembl ids as the identifier? switch to this??

# make sure raw counts are in spatial X slot 
if "counts" in ST_adata.layers:
    ST_adata.X = ST_adata.layers["counts"]
else:
    ST_adata.X = ST_adata.X 

# remove mitochondria-encoded (MT) genes from spatial if present
# mito gene symbols start with "MT-"
ST_adata.var['MT_gene'] = [gene.startswith('MT-') for gene in ST_adata.var['SYMBOL']]
ST_adata.obsm['MT'] = ST_adata[:, ST_adata.var['MT_gene'].values].X.toarray()
ST_adata = ST_adata[:, ~ST_adata.var['MT_gene'].values]

# check for duplicated gene ids
dup = ST_adata.var.index.duplicated(keep='first')
ST_adata.var['dup'] = dup
ST_adata = ST_adata[:, ~ST_adata.var['dup'].values]

print("Training single cell model")
cell2location.models.RegressionModel.setup_anndata(
    adata=adata_sc,
    batch_key=sc_batch,
    labels_key=sc_annotation_column,
    categorical_covariate_keys=sc_categorical_covariate
    )

# create the regression model
mod = cell2location.models.RegressionModel(adata_sc)
mod.train(
    max_epochs=args.epochs_ref,
    batch_size=args.sc_batch_size,
    accelerator="gpu"
    )
print("GPU usage after training single cell model")
print_gpu_usage()

print("Single cell loss figure")
plt.figure()
mod.plot_history(20)
plt.legend(labels=['scRNASeq']);
plt.savefig(f"{reports_folder}/loss_scRNASeq.png")

adata_sc = mod.export_posterior(
    adata_sc,
    sample_kwargs={'num_samples': 1000, 'batch_size': args.sc_batch_size, 'accelerator': 'gpu'}
)

print("GPU usage after exporting single cell posterior")
print_gpu_usage()

# Save model
mod.save(f"{ref_model}", overwrite=True)

print("RNA-seq QC plot")
plt.figure()
mod.plot_QC()
plt.title("QC scRNAseq")
plt.savefig(f"{reports_folder}/QC_scRNASeq.png")

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_sc.varm.keys():
    inf_aver = adata_sc.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_sc.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_sc.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_sc.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_sc.uns['mod']['factor_names']
print(f"\nEstimated expression per cluster{inf_aver.iloc[0:2,:]}\n")

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(ST_adata.var_names, inf_aver.index)
print('\nIntersect with sc:', len(intersect), "genes\n")

if len(intersect) == 0:
    inf_aver.index = adata_sc.var.feature_name.values
    intersect = np.intersect1d(ST_adata.var_names, inf_aver.index)
    print('\nChanging to symbol Intersect with sc:', len(intersect), "genes\n")
    
ST_adata = ST_adata[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

print(type(inf_aver))
inf_aver.to_csv(os.path.join(output_folder, "reference_cell_type_signature.csv"), index = True)

print("GPU usage before training spatial model")
print_gpu_usage()

print(torch.cuda.memory_summary(device=None, abbreviated=False))

torch.cuda.empty_cache()
print(torch.cuda.memory_summary(device=None, abbreviated=False))

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(
    adata=ST_adata,
    batch_key=st_batch,
    categorical_covariate_keys=st_categorical_covariate
    )

mod = cell2location.models.Cell2location(
    ST_adata,
    cell_state_df=inf_aver, 
    N_cells_per_location=8,
    detection_alpha=20
    )

print(adata_sc)
print(ST_adata)

print("\nSTART TRAINING\n")
mod.train(max_epochs=args.epochs_st,
          batch_size=args.st_batch_size,
          train_size=1,
          accelerator="gpu"
         )

print("\nTraining done\n")
print(mod.adata.n_obs)

print("plot ELBO loss history during training, removing first 100 epochs from the plot")
plt.figure()
mod.plot_history(100)
plt.legend(labels=['full ST data training']);
plt.savefig(f"{reports_folder}/loss_ST.png")

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
if args.st_batch_size == None:
    batch_size = mod.adata.n_obs
else:
    batch_size = args.st_batch_size
ST_adata = mod.export_posterior(
    ST_adata,
    sample_kwargs={'num_samples': 1000, 'batch_size': batch_size, 'accelerator': 'gpu'}
    )

# Save model
mod.save(f"{st_model}", overwrite=True)

plt.figure()
mod.plot_QC()
plt.title("Quality control ST")
plt.savefig(f"{reports_folder}/QC_ST.png")


print("\ncreating csv\n")
# add 5% quantile, representing confident cell abundance, 'at least this amount is present', to adata.obs with nice names for plotting
results = pd.DataFrame(ST_adata.obsm['q05_cell_abundance_w_sf'])
results.columns = ST_adata.uns['mod']['factor_names']

# normalize rows from 0 to 1, convert to percentage
results = results.div(results.sum(axis=1), axis=0)
results.to_csv(os.path.join(output_folder, 'cell_abundance_q5.csv'), index=True)
