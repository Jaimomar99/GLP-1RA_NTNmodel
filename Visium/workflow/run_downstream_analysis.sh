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


Rscript ${PIPELINE}/R/downstream_analysis.R \
    --manifest_path ${SAMPLE_MANIFEST} \
    --deconvolute_out ${DECONVOLUTION_OUT} \
    --deconvolute_nsplits ${DECONVOLUTE_NSPLITS} \
    --integrated_out ${INTEGRATED_OUT} \
    --by_sample_out ${BY_SAMPLE_OUT} \
    --ultimate_processed_out ${ULTIMATE_PROCESSED_OUT} \
    --loupe_browser_out ${LOUPE_BROWSER_OUT} \
    --pipeline_path ${PIPELINE} \
    --feature_df_path ${FEATURE_METADATA} \
    > ${BASE_DIR}/${DATASET_NAME}/workflow/logs/downstream_analysis${DATE}.log 2>&1
 