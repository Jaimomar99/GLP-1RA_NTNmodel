#!/bin/bash

# Parse command line arguments
while getopts ":c:i:" opt; do
  case $opt in
    c)
      config_file=$OPTARG
      ;;
    i)
      array_id=$OPTARG
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

# change directory to where output will be written
cd ${PROJECT_DIR}/processed/spaceranger_blank


IFS=','; read -r ALIAS_ID SEQ_SAMPLE_ID FW_SAMPLE_ID VISIUM_SLIDE VISIUM_CAPTURE_AREA DUAL_INDEX ELN DATE_PROCESSED_LAB DATE_SEQUENCING SEQUENCING_RUN DIAGNOSIS DONOR_DESCRIPTION IMAGE <<< "$(sed -n ${array_id}p "${SAMPLE_MANIFEST}")"
    echo "Running spaceranger on sample $ALIAS_ID"
    ${SPACERANGER} count \
        --id=${ALIAS_ID}_blank \
        --sample=${SEQ_SAMPLE_ID} \
        --description=${DIAGNOSIS} \
        --transcriptome="${TRANSCRIPTOME}" \
        --probe-set="${PROBESET}" \
        --fastqs="${PROJECT_DIR}/raw/fastq/${SEQUENCING_RUN}" \
        --image="${BLANK_VISIUM_IMAGE}" \
        --slide=${VISIUM_SLIDE} \
        --area=${VISIUM_CAPTURE_AREA} \
        --jobmode local \
        --localcores=32 \
        --localmem=128
    echo "sample ${ALIAS_ID} finished running"
