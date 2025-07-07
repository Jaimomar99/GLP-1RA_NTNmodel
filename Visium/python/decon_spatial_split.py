#!/usr/bin/env python

import scanpy as sc
import pandas as pd
import os
import argparse
from sklearn.model_selection import train_test_split


def stratified_split(df, target_col, num_splits):
    splits = []
    split_size = 1.0 / num_splits

    for _ in range(num_splits - 1):
        print(split_size)
        df, split = train_test_split(df, test_size=split_size, stratify=df[target_col])
        splits.append(split)
        # Update split_size and df for the next iteration
        split_size = split_size / (1.0 - split_size)

    splits.append(df)

    return splits


parser = argparse.ArgumentParser()

parser.add_argument("--project_folder", type=str, help="Indicate the path to the project")
parser.add_argument("--out_dir", type=str, help="output_subdirectory in processed")
parser.add_argument("--st_path", type=str, help="Indicate the path with the spatial h5ad object")
parser.add_argument('--deconvolute_split', type=bool, default=False, help='Whether to run the split function. Should be "True" or "False".')
parser.add_argument('--target_cols', type=str, help='Column to stratify by.')
parser.add_argument('--num_splits', type=int, help='Number of splits.')

args = parser.parse_args()

# SET VARIABLES
output_folder = f"{args.project_folder}/processed/{args.out_dir}/"

if args.deconvolute_split:
    st_adata = sc.read_h5ad(args.st_path)

    # this is necessary only when using Seurat to write h5ad, move to integration script that write combined h5ad?
    st_adata._raw._var.rename(columns={'_index': 'features'}, inplace=True)

    # which columns to use to stratify by in split
    stratify_list = [x.strip() for x in args.target_cols.split(",")]
    print(stratify_list)

    spot_df = st_adata.obs.copy()
    spot_df["strat_cat"] = spot_df[stratify_list].apply(tuple, axis=1)

    split_list = stratified_split(spot_df, target_col="strat_cat", num_splits=args.num_splits)

    os.makedirs(f"{output_folder}/split_h5ad", exist_ok=True)
    for i, s in enumerate(split_list):
        split_adata = st_adata[s.index, ]
        split_adata.write_h5ad(filename=f"{output_folder}/split_h5ad/st_split_{i}.h5ad")
