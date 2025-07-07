#!/usr/bin/env python

import torch
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import cell2location
from cell2location.utils.filtering import filter_genes
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

if len(args.sc_covariates) == 0:
    sc_categorical_covariate = None
else:
    covariate_list = args.sc_covariates.split(",")
    sc_categorical_covariate = [x.strip() for x in covariate_list]

adata_sc = sc.read_h5ad(sc_path)

# finish setting up single cell and spatial anndata objects
adata_sc.var['SYMBOL'] = adata_sc.var_names # feature ids are gene symbols

# make sure raw counts are in single cell X slot 
if "counts" in adata_sc.layers:
    adata_sc.X = adata_sc.layers["counts"]
else:
    # assume adata_sc.X is the raw counts
    pass
 
# filter genes
selected = filter_genes(adata_sc, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_sc = adata_sc[:, selected].copy()


print("Training single cell model")
print(sc_batch)
print(sc_annotation_column)
print(sc_categorical_covariate)
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

# Save anndata object with results
del adata_sc.raw  # to avoid errors when writing h5ad (happens w/ h5ad created w/ seuratdisk)
adata_file = f"{ref_model}/sc.h5ad"
adata_sc.write(adata_file)

print("RNA-seq QC plot")
plt.figure()
mod.plot_QC()  # plots are saved on top of each other in png; open issue in: https://github.com/BayraktarLab/cell2location/issues/191
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

print("Writing reference cell type signature to file")
inf_aver.to_csv(os.path.join(ref_model, "reference_cell_type_signature.csv"), index = True)
