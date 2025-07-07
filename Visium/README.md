# visium_pipeline

Analysis pipeline for processing 10X Visium data sets.

## Table of Contents

- [Environment setup](#environment-setup)
- [Cell ranger processing](#cell-ranger-processing)
- [QC](#qc)
- [workflow](#workflow)
- [License](#license)

## Environment setup 


### Prerequisites

- Spaceranger version 2.1.1:

`wget -O spaceranger-2.1.1.tar.gz "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-2.1.1.tar.gz?Expires=1696397354&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=JDMZcQBQokLsdN4H~x-Dwll7tTUIaa-31lAY1Mgf7UQh5cJ1dmxssZF10-0Xy0P7b5Tzpm4byw2oPPFje2HRccIhqSc9cfWM5L2pAPL-oLMPvs7JsDDknmO~QcBgAiwDZlrQJGnAC6qZl2~tYB-9Tf5wk0UEf8MLc47Bg5hYs7ScIm2mho7ABdo7oq3eS0Gwn29KuRjJ42rPDjstbsus3co-xpjKtd33vv4V7AGSAAI4cncNadkb8KUzNLSfRdtYCyAoJoP8qi4NDw-XRjIFH2~RGaJ0FsiDTm~Nhu9W8cNmO9UnOTMT3Twd53thQQPjUj1y1vjv4p-PHxvxCHicbQ__"`

- Spaceranger GRCh38 reference:

`wget "https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-GRCh38-2020-A.tar.gz"`

- Spaceranger mm10 reference:

`wget "https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-mm10-2020-A.tar.gz"`


## workflow
|          |        Workflow               |                                definition                                  |
|----------|-------------------------------|----------------------------------------------------------------------------|
| Platform | 0. sh Initate project         | creating folders                                                           |
|          | 1. sh Space ranger            | from fastq to matrices                                                     |
|          | 2. Rscript by sample          | QC, SCT, umap, clustering..                                                |
|          | 3. Integration & DE           | merge,  harmony and DE by cluster                                          |
|          | 4. python deconvolution       | cell2location                                                              |
|          | 6. Mapping & generation       | map clustering and abundance to slides  & create counts, metadata and umap |

### Initiate project
1. Create manifest csv file with all samples that are part of dataset; see `workflow/templates/manifest.csv` for an example template.
Internally generated datasets that are processed from raw sequencing data should include the columns listed below. Public datasets may not need all of the columns (ex. if spaceranger output is donwloaded as starting point)
   - alias_id : unique identifier (integer or character string) to use in place of FW sample id 
   - seq_sample_id : this is in the sequencing sample worksheet and is necessary for demultiplexing
   - fw_sample_id : FW sample identifier; may not be applicable for external/public datasets or non-human datasets
   - visium_slide : Visium slide id (ex. V12N14-315); this is used when running spaceranger
   - visium_capture_area : capture are on Visium slide (ex. B1)
   - dual_index : sequencing dual index; this is in the sequencing sample worksheet
   - date_processed_lab : date experiment was initiated in the lab (ex. Visium slides and library prep)
   - date_sequencing : date of sequencing run
   - sequencing_run : sequencing run id (ex. HGWLLDMXY); this is necessary for demultiplexing/locating fastq
   - diagnosis : disease state or other categorical grouping variable (ex. T2D vs. healthy)
   - donor_description : donor description
   - image : absolute path to high resolution image that is used as input to spaceranger
   - manual_align_json: if manual alignment to fiducial markers is needed; store path to json file here

   ADD: version (ex. FF, FFPE)
   
2. Make dataset.config; see template in `workflow/templates/manifest.csv`
3. Activate the correct shared conda environment b/f running any of the pipeline steps
4. Set up the dataset directory by running: `./setup_dataset.sh -c dataset.config`
   - Move dataset.config into project workflow folder
5. Stage raw data for analysis using custom `stage_data.sh` script:
   - Copy high resolution images to `raw/images`; these images are used as input to spaceranger
   - Copy fastq files to `raw/fastq`; fastq should be organized by sequencing library and/or sequencing run if there was more than one
      - only demultiplexed sequencing data should be stored in `raw/fastq`; demultiplexing is not part of this pipeline, so demultiplex separately and stage only fastq files

### Run Spaceranger
6. Run spaceranger to generate filtered count matrix and spot spatial coordinates:
   - `workflow/templates/spaceranger_slurm.sh` is a template job script for submitting spaceranger array jobs; it runs the `run_spaceranger.sh` script
   - Copy `spaceranger_slurm.sh` template into `workflow` and modify as needed (ex. change dataset.config and set SLURM_ARRAY_TASK_ID to indicate array job indeces)
      - array job indeces refer to the sample row in the manifest csv
   - `sbatch spaceranger_slurm.sh` => this submits an array job to slurm; one job for each sample in the manifest
7. Go over spaceranger output:
   - Use open on demand desktop to look at the html summary from spaceranger
   - Locally: copy html summary files to your local laptop/desktop
   - Do rough sample level QC check; (ex. # of genes per spot; total # of UMIs, UMIs under tissue)
   - Check fiducial marker alignment and note which samples need to be manually aligned
      - manual alignment can only be done locally in Loupe browser; see manual alignment documentation here:
      - download high resolution image to local laptop/desktop and note Visium slide and capture area before doing alignment workflow in Loupe browser
      - after completing manual alignment, make sure to save json file, upload back to marjorie, and add absolute path to manifest csv in `manual_align_json` column
   - Delete initial spaceranger run, and re-run spaceranger for any samples that needed manual fiducial alignment

### Analysis
8. Run by sample analysis (run_analysis.sh) : create slurm script
9. Run integration (now also in run_analysis.sh): create slurm script
10. Run differential expression ()
   - Wilcoxon rank sum: comparing clusters (automated as part of pipeline)
   - Comparing disease state (ex. healthy vs disease per cluster/cell type)
11. Curate single cell datasets for use as deconvolution reference (ex. in jupyter notebook)
   - Check that state/condition matches the spatial data you want to deconvolute
   - Check basic QC, 
   - Save subset-ed anndata as h5ad
     - Make note of SC_BATCH_KEY or SC_COVARIATES (correction for decon)
     - Add to documentation (ex. README)
12. Generate deconvolution reference
   - Copy slurm template to project workflow folder (deconvolution_slurm_ref.sh)
   - Run_deconvolution_ref.sh using template in deconvolution_slurm_ref.sh
     - Sbatch decon…
   - Check QC plots for model performance
13. Run deconvolution using reference
   - Split spatial samples if necessary (ex. if > 24,000 spots, run split deconvolution)
     - Modify dataset_config (number of splits)
   - Copy slurm template to project workflow folder (deconvolution_slurm.sh)
   - Run_deconvolution.sh from template
14. Run downstream analysis (run_downstream_analysis.sh)
   - Map deconvolution results to integrated object
   - Calculate Compositional cluster
   - Map integration clustering, compositional cluster and deconv to each sample
   - Create loupe browsers
   - Cell type plots to check proportions, either invididually or pie chart
