#!/bin/bash

# Parse command line arguments
while getopts ":c:" opt; do
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

# Read variables from config file

source "${config_file}"
source "${CONDA_PROFILE}"
conda activate "${CONDA_ENV}"

# setup to read variables from config file
DATE=`date +"%Y-%m-%d_%b"`
PROJECT_DIR="${BASE_DIR}/${DATASET_NAME}"

echo $BASE_DIR

# TO DO: add some status/progress messages
# set up project directory
Rscript ${PIPELINE}/R/0_setup_dataset.R \
    --dataset_name ${DATASET_NAME} \
    --dir_path ${BASE_DIR}

# To do: check if they exist; don't want to overwrite
cp "${PIPELINE}/workflow/templates/readme_template.md" ${PROJECT_DIR}/README.md
cp "${PIPELINE}/workflow/templates/spaceranger_slurm.sh" ${PROJECT_DIR}/workflow