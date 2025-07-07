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

PROJECT_DIR="${BASE_DIR}/${DATASET_NAME}/"
PROCESSED_DIR="${BASE_DIR}/${DATASET_NAME}/processed/"

### by sample notebook
PLOT_DIR=${PROJECT_DIR}/notebooks/${BY_SAMPLE_OUT}/plots/
OUT_DIR=${PROJECT_DIR}/notebooks/${BY_SAMPLE_OUT}

papermill \
    ${PIPELINE}/workflow/templates/notebooks/by_sample.ipynb \
    ${PROJECT_DIR}/notebooks/${BY_SAMPLE_OUT}.ipynb \
    -p plot_dir "${PLOT_DIR}" \
    -p outdir "${OUT_DIR}" \
    -p project_dir "${PROJECT_DIR}" \
    -p pipeline_dir "${PIPELINE}" \
    -p processed_dir "${PROCESSED_DIR}" \
    -p by_sample_out "${PROCESSED_DIR}/${BY_SAMPLE_OUT}/" \
    -p manifest_path "${SAMPLE_MANIFEST}" \
    -p write_plots FALSE \
    --prepare-only

### integrated notebook
PLOT_DIR=${PROJECT_DIR}/notebooks/${INTEGRATED_OUT}/plots/
OUT_DIR=${PROJECT_DIR}/notebooks/${INTEGRATED_OUT}

papermill \
    ${PIPELINE}/workflow/templates/notebooks/integrated.ipynb \
    ${PROJECT_DIR}/notebooks/${INTEGRATED_OUT}.ipynb \
    -p plot_dir "${PLOT_DIR}" \
    -p outdir "${OUT_DIR}" \
    -p project_dir "${PROJECT_DIR}" \
    -p pipeline_dir "${PIPELINE}" \
    -p manifest_path "${SAMPLE_MANIFEST}" \
    -p integrated_rds "${PROCESSED_DIR}/${INTEGRATED_OUT}/harmonized.rds" \
    -p write_plots FALSE \
    --prepare-only
