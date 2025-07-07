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

# Reference celltype signatures
${PIPELINE}/python/cell2location_reference.py \
  --project_folder ${BASE_DIR}/${DATASET_NAME} \
  --out_dir "${DECONVOLUTION_OUT}" \
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
  > ${BASE_DIR}/${DATASET_NAME}/workflow/logs/deconvolution_ref_${DATE}.log 2>&1
