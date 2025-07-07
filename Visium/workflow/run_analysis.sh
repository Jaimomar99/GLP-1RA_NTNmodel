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


Rscript ${PIPELINE}/R/2_by_sample.R \
    --manifest_path ${SAMPLE_MANIFEST} \
    --out_dir ${BY_SAMPLE_OUT} \
    --pipeline_path ${PIPELINE}\
    > ${BASE_DIR}/${DATASET_NAME}/workflow/logs/by_sample_${DATE}.log 2>&1


Rscript ${PIPELINE}/R/3_integration.R \
    --manifest_path ${SAMPLE_MANIFEST} \
    --by_sample_subdir ${BY_SAMPLE_OUT} \
    --out_dir ${INTEGRATED_OUT} \
    --grouping ${GROUP_VAR} \
    --batch_correction_var "${BATCH_CORRECTION}" \
    --pipeline_path ${PIPELINE}\
    > ${BASE_DIR}/${DATASET_NAME}/workflow/logs/integration_${DATE}.log 2>&1
 

Rscript ${PIPELINE}/R/DE_clustering.R \
    --manifest_path ${SAMPLE_MANIFEST} \
    --integrated_out ${INTEGRATED_OUT} \
    --pipeline_path ${PIPELINE}\
    > ${BASE_DIR}/${DATASET_NAME}/workflow/logs/de_clustering_${DATE}.log 2>&1