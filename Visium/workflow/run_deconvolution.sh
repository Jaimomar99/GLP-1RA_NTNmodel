#!/bin/bash

# Parse command line arguments
while getopts ":c:i:" opt; do
  case $opt in
    c)
      config_file=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Check if config file is provided
if [ -z "$config_file" ]; then
  echo "Please provide a path to the config file using the -c option."
  exit 1
fi

# load cuda to use gpu
module load cuda/11.8.0-spk

# Read variables from config file
source "${config_file}"
source "${CONDA_PROFILE}"
conda activate ${CONDA_ENV}
DATE=`date +"%Y-%m-%d_%b"`


if [ "$DECONVOLUTE_SPLIT" = "True" ]; then
    # split spatial data for deconvolution if necessary
    ${PIPELINE}/python/decon_spatial_split.py \
        --project_folder ${BASE_DIR}/${DATASET_NAME} \
        --out_dir "${DECONVOLUTION_OUT}" \
        --st_path ${ST_PATH} \
        --deconvolute_split ${DECONVOLUTE_SPLIT} \
        --target_cols "${DECONVOLUTE_SPLIT_STRATIFY}" \
        --num_splits ${DECONVOLUTE_NSPLITS} \
        > ${BASE_DIR}/${DATASET_NAME}/workflow/logs/split_for_decon_${DATE}.log 2>&1

    # run deconvolution on each split of the data
    echo "running on split h5ad's"
    NSPLITS=$((DECONVOLUTE_NSPLITS - 1))
    for i in $(seq 0 ${NSPLITS})
    do
        CURRENT_ST_PATH="${BASE_DIR}/${DATASET_NAME}/processed/${DECONVOLUTION_OUT}/split_h5ad/st_split_${i}.h5ad"

        ${PIPELINE}/python/cell2location_decon_spatial.py \
            --project_folder ${BASE_DIR}/${DATASET_NAME} \
            --out_dir "${DECONVOLUTION_OUT}/split_${i}" \
            --reference_celltypes ${BASE_DIR}/${DATASET_NAME}/processed/${DECONVOLUTION_OUT}/reference_signatures/reference_cell_type_signature.csv \
            --st_path ${CURRENT_ST_PATH} \
            --sc_path ${SC_PATH} \
            --sc_batch ${SC_BATCH_KEY} \
            --st_batch ${ST_BATCH_KEY} \
            --sc_covariates "${SC_COVARIATES}" \
            --st_covariates "${ST_COVARIATES}" \
            --sc_annotation_column ${SC_ANNOTATION} \
            --sc_batch_size ${SC_BATCH_SIZE} \
            --st_batch_size ${ST_BATCH_SIZE} \
            --epochs_ref ${TRAIN_EPOCHS_REF} \
            --epochs_st ${TRAIN_EPOCHS_ST} \
            --cells_per_spot ${CELLS_PER_SPOT} \
            > ${BASE_DIR}/${DATASET_NAME}/workflow/logs/deconvolution_split_${i}_${DATE}.log 2>&1
    done
else
    echo "running on full h5ad without splitting"
    ${PIPELINE}/python/cell2location_decon_spatial.py \
    --project_folder ${BASE_DIR}/${DATASET_NAME} \
    --out_dir "${DECONVOLUTION_OUT}" \
    --reference_celltypes ${BASE_DIR}/${DATASET_NAME}/processed/${DECONVOLUTION_OUT}/reference_signatures/reference_cell_type_signature.csv \
    --st_path ${ST_PATH} \
    --sc_path ${SC_PATH} \
    --sc_batch ${SC_BATCH_KEY} \
    --st_batch ${ST_BATCH_KEY} \
    --sc_covariates "${SC_COVARIATES}" \
    --st_covariates "${ST_COVARIATES}" \
    --sc_annotation_column ${SC_ANNOTATION} \
    --sc_batch_size ${SC_BATCH_SIZE} \
    --st_batch_size ${ST_BATCH_SIZE} \
    --epochs_ref ${TRAIN_EPOCHS_REF} \
    --epochs_st ${TRAIN_EPOCHS_ST} \
    --cells_per_spot ${CELLS_PER_SPOT} \
    > ${BASE_DIR}/${DATASET_NAME}/workflow/logs/deconvolution_${DATE}.log 2>&1
fi
