This repo contains all code to reproduce analysis and figures in "Mapping cell-cell fusion at single-cell resolution" (https://doi.org/10.1101/2024.12.11.627873)

### Accessing raw data
Raw and processed single cell RNA-sequencing can be downloaded from GEO with accession number GSE286213. Targeted barcode sequencing to generate a list of barcodes in the samples using pycashier can be downloaded from Gene Expression Omnibus (GEO) with accession number GSE286212.

In notebooks 7-13, the `h1806_colors_cellclass_hclust_nb-barassigned_ccclustered_precc_singlets.rds` is equivalent to 
`GSE286213_CMDuo_seurat_qc_barcodeAnnotated_HCC1806_samples.rds` from GEO, and `mb231_colors_cellclass_hclust_nb-barassigned_ccclustered_precc_singlets.rds` is equivalent to `GSE286213_CMDuo_seurat_qc_barcodeAnnotated_MDAMB231_samples.rds`. Start here if you want to perform downstream analysis. These objects can also be used as the starting point in notebooks 3-6, the objects in notebooks 1-2 need to be built from raw. 

## Raw data processing
To start from raw sequencing data, follow the notes in `00_Raw_data_processing_and_alignment_notes.md`


### Setting up the Rstudio analysis environment
To setup your analysis environment, download the Docker Image which is equipped to launch Rstudio with all the necessary packages and versions:
`docker pull phdidi/cmduo_rocker:latest`

Port Rstudio to your local browser from your server by running:
`docker run --rm -ti -e PASSWORD=makeapassword -e USERID="$(id -u)" -e GROUPID="$(id -g)" -p 8770:8787 -v "$(pwd)":/home/rstudio/workspace cmduo_rocker`

On your local machine run (change 'user' to your username on your server and 'host' to your server name)
`ssh -NL 8770:localhost:8770 -p 22 user@host`

Then open a browser and navigate to: http://localhost:8770/⁠

When prompted, enter 'rstudio' as your username and 'makeapassword' (or whatever you set this to) as your password to launch Rstudio.



