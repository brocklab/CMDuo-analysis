Start at this notebook to reproduce analysis from raw data. Many of the steps in this notebook are computationally-intensive and will take time. Unless there is a need to perform analysis from the raw data (i.e. mapping to a different reference genome), it is recommended that you instead start with analysis at the stage of Seurat creation.

Unless otherwise stated, these steps are performed from the command line in a Linux environment.

This notebook walks through a few tasks:
1. Data is imported from four different sources
	1. Parse whole transcriptome Illumina sequencing
	2. Parse CRISPR detect Illumina sequencing
	3. Targeted barcode Illumina sequencing
	4. Reference genome and annotations
2. Pycashier is run on targeted barcode sequencing to generate a list of ClonMapper barcodes present in the samples
3. Pycashier outputs are processed in R to create a Parse CRISPR Detect pipeline compatible list of guides/barcodes
4. The Parse pipeline is installed
5. GFP and mCherry are added to the genome and a reference genome is built
6. Single-cell whole transcriptome data is aligned to the reference genome for each sublibrary then combined
7. Single-cell CRISPR detect data is aligned to the pycashier list of barcodes for each sublibrary then combined

## 1. Download and organize data (update when GEO is ready)

#### i. Raw data from Parse WT and CRISPR Detect
This is the data generated from scRNA-seq using the Parse WT kit (v2 chemistry, cat. no. ECW02130) with the CRISPR-detect add-on (cat. no. CRS1010). For each sample, Read 1 and read 2 fastq files from multiple lines have been pre-concatenated into the single file stored on GEO for each sample and read. 
Note: Sample names are 1-4 and 6-8. A pipetting error occurred during early laboratory preparation of Parse sublibrary 5, so it was excluded from further preparation and analysis.

i. Download the raw Parse whole transcriptome fastqs from GEO and place in a folder called WT. 
```
mkdir expdata

#to-do
```

in 'parsealign/expdata/WT', these files should be present for each sample (1-4, 6-8): 
WT-{sample_num}.R1.raw.fastq.gz  
WT-{sample_num}.R2.raw.fastq.gz

#### ii. Download the raw Parse CRISPR detect fastqs from GEO and place in a folder called CRISPRdetect
```
#to-do
```

in 'parsealign/expdata/WT', these files should be present for each sample (1-4, 6-8): 
CR-{sample_num}.R1.raw.fastq.gz  
CR-{sample_num}.R2.raw.fastq.gz

#### iii. Raw data from targeted barcode sequencing
The ClonMapper and ClonMapper Duo technologies express genomic barcodes as sgRNAs and polyadenylated transcripts, allowing detection of barcodes through multiple means. The Parse CRISPR-detect pipeline was designed to work with CROPseq guides for CRISPR screens, but is easily adaptable to the ClonMapper platform. In order to use CRISPR detect with ClonMapper, we will need to provide the Parse pipeline a list of guides/barcodes to align to. As ClonMapper barcodes are randomly generated, targeted Illumina sequencing (https://docs.brocklab.com/clonmapper/protocol/#clonmapper-barcode-sampling-of-cells) is be performed to identify the set of barcodes contained within each sample. 

i. Download targeted sequencing data from GEO
```
mkdir expdata/targeted_barcode

to-do
```

```
seqdat move JA24019 --out ./expdata/targeted_barcode
```
#### iv. Reference genome and annotations
Transcriptomic sequencing data always needs to be aligned to a reference genome. To generate reference annotations, it is recommended that users use to the latest version of the genomic annotations available for their organism. Since we have used human cancer cell lines in this study, we will use the human reference genome, GRCh38. 

a. Get the latest version of the human genome (sequence and annotation files) from Ensembl (at the time this notebook was generated, the current release was v112). Note: This process may take up 20 minutes or more depending on the speed of your internet connection
```
mkdir -p parsealign/genomes
cd parsealign

wget -P genomes https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz;
wget -P genomes https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz 
```

## 2. Run pycashier to generate a list of ClonMapper barcodes
The ClonMapper and ClonMapper Duo technologies express genomic barcodes as sgRNAs and polyadenylated transcripts, allowing detection of barcodes through multiple means. The Parse CRISPR-detect pipeline was designed to work with CROPseq guides for CRISPR screens, but is easily adaptable to the ClonMapper platform. In order to use CRISPR detect with ClonMapper, we will need to provide the Parse pipeline a list of guides/barcodes to align to. As ClonMapper barcodes are randomly generated, targeted Illumina sequencing (https://docs.brocklab.com/clonmapper/protocol/#clonmapper-barcode-sampling-of-cells) is performed to identify the set of barcodes contained within each sample.

Note:  The Parse CRISPR Detect pipeline will not run if you provide if a list of unique barcodes in which any pair of barcodes has a Hamming distance < 3. Barcode pairs will need to be collapsed.

#### i. Combine fastqs into early time point (initial and pre-sort) and late time point (fusion and control) 
Doing this so that starcode clusters across many samples and so that samples can be sorted based on initial abundance in the early timepoint samples
```
mkdir pycash_analysis_clean

cat *initial*.fastq.gz *presort*.fastq.gz > combined/combined_H1806_early.fastq.gz
cat *fusion*.fastq.gz *ctrl*.fastq.gz > combined/combined_H1806_late.fastq.gz
```

#### ii. Run pycashier on targeted sequencing data. 
A docker implementation of pycashier (https://github.com/brocklab/pycashier) can be used to extract barcodes from targeted sequencing data. Here, we updated the distance to 2 and the clustering ratio to 1 in order to cluster similar barcodes.
```
docker run --rm -it -v $PWD:/data -u $(id -u):$(id -g) ghcr.io/brocklab/pycashier:v2024.1005 extract \
-i ./expdata/targeted_barcode/h1806/combined \
-o ./pycash_analysis_clean/pycash_outs/h1806 \
-d 2 \
-r 1 \
--threads 16

docker run --rm -it -v $PWD:/data -u $(id -u):$(id -g) ghcr.io/brocklab/pycashier:v2024.1005 extract \
-i ./expdata/targeted_barcode/mb231/combined \
-o ./pycash_analysis_clean/pycash_outs/mb231 \
-d 2 \
-r 1 \
--threads 16
```

#### iii. Combine barcodes detected across all samples into one file for each cell line
Running pycashier receipt combines multiple pycashier outs and perserves sample information in a new column
```
docker run --rm -it -v $PWD:/data -u $(id -u):$(id -g) ghcr.io/brocklab/pycashier:v2024.1005 receipt \
-i ./pycash_analysis_clean/pycash_outs/h1806 \
-o ./pycash_analysis_clean/combined_d2r1_h1806_only.tsv

docker run --rm -it -v $PWD:/data -u $(id -u):$(id -g) ghcr.io/brocklab/pycashier:v2024.1005 receipt \
-i ./pycash_analysis_clean/pycash_outs/mb231 \
-o ./pycash_analysis_clean/combined_d2r1_mb231_only.tsv
```

## 3. Process barcodes to format input list for Parse CRISPR Detect pipeline

Here, we perform post-processing of the barcodes in R to remove noise, double check barcode similarity, and name barcodes based on their color and overall abundance in the early timepoint samples. Once this is done, we format the named barcodes from each cell line into a single file for use during the Parse CRISPR detect alignment step.

```{r}
library(tidyverse)
library(stringdist)
`%notin%` <- negate(`%in%`)

## read in data
pycashbars_h1806 <- read.csv('combined_d2r1_h1806_only.tsv', sep='\t', header=TRUE)
pycashbars_mb231 <- read.csv('combined_d2r1_mb231_only.tsv', sep='\t', header=TRUE)

# barcodes in both samples are likely to be noise or contamination. Detect these for filtering
shared <- intersect(pycashbars_h1806$barcode, pycashbars_mb231$barcode)

format_barcodes_cr_detect <- function(pycash_outs, cellline, min_lev = 3){
  # barcodes in both samples are likely to be noise or contamination, remove these
  pycash_outs <- pycash_outs[pycash_outs$barcode %notin% shared,]

  # set GFP and mCherry index sequences
  gfp_indices <- c('GACAA','GCTGT','ATCGC','ACGCA')
  mch_indices <- c('CATCC','CGAGA','TAGTG','TCAAC')
  
  # assign barcode color to each detected guide, by matching to index within a Levenshtein distance of 1
  named_guides <- pycash_outs %>% 
    separate(sample, into = c("datatype", "cellline", "condition"), sep = "_") %>% 
    group_by(barcode, condition) %>% 
    summarize(condition_count = sum(count)) %>% 
    arrange(desc(condition_count)) %>% 
    arrange(condition) %>% select(barcode) %>% 
    distinct() %>% ungroup() %>% 
    mutate(barnum = 1:n())  %>% 
    mutate(ind=substr(barcode, 1, 5)) %>% 
    rowwise() %>% 
    mutate(color = case_when(any(stringdist(ind, c(gfp_indices)) <= 1) & any(stringdist(ind, c(mch_indices)) > 1)  ~ "gfp",
                             any(stringdist(ind, c(mch_indices)) <= 1) & any(stringdist(ind, c(gfp_indices)) > 1) ~ "mcherry",
                             .default="nocolor")) %>% 
    mutate(barid = case_when(color=='gfp'     ~ paste0('g', str_pad(barnum, 3, pad = "0")),
                             color=='mcherry' ~ paste0('r', str_pad(barnum, 3, pad = "0")),
                             color=='nocolor' ~ paste0('u', str_pad(barnum, 3, pad = "0")),
                             .default = NA) ) %>% ungroup() 
  
  # check to see if any barcodes are highly similar via Levenshtein distance
  lv_mat <- stringdistmatrix(named_guides$barcode, named_guides$barcode, method = "lv")
  similar_bars <- length(lv_mat[lv_mat <= min_lev][lv_mat[lv_mat <= min_lev] != 0])/2
  
  print(paste(length(named_guides$barcode), 'barcodes detected in', cellline,'cells.', sep=' '))
  
  if(similar_bars > 0){
    print(paste('WARNING: there are', similar_bars,'pairs of barcodes with Levenshtein distance',min_lev,'or less.', sep=' '))
  }
  
  if(similar_bars == 0){
    print(paste('Looks good: there are no pairs of barcodes with Levenshtein distance',min_lev,'or less.', sep=' '))
  }
  
  # Final formatting for Parse
  parse_guides_pycashbars <- data.frame(
        Guide_Name=paste(named_guides$barid, cellline, sep='_'), 
    Guide_Sequence=named_guides$barcode, 
       Target_Gene=named_guides$barcode) %>% 
    mutate(Prefix='ATCTTGTGGAAAGGACGAAACACCG',
           Suffix='GTTTTAGAGCTAGAAATAGCAAGTT') %>% 
    select(Guide_Name, Prefix, Guide_Sequence, Suffix, Target_Gene)
  
  return(parse_guides_pycashbars)
}

formatted_h1806 <- format_barcodes_cr_detect(pycashbars_h1806, cellline='h1806', min_lev=2)
formatted_mb231 <- format_barcodes_cr_detect(pycashbars_mb231, cellline='mb231', min_lev=2)


final_out <- rbind(formatted_h1806, formatted_mb231)
write.table(final_out, 'formatted_guide_list_cellline_named.csv', sep=",", row.names = FALSE, quote=FALSE)
```


## 4. Setup environment and install Parse pipeline

Most of this code and workflow is standard and adapted from Parse
#### i. Create and activate a conda environment with Python version >=3.1 
```
conda create -n spipe conda-forge::python==3.10  
conda activate spipe
```
#### ii. Install Parse pipeline (v 1.3.1)
This needs to be downloaded from the Parse Support Suite. Can scp the local download to your server analysis folder.
```
scp -P 22 ParseBiosciences-Pipeline.1.3.1.zip {username}@{server-address}:{path-to-folder}/parsealign
```
#### iii. Unzip the file on your server, navigate into unzipped folder
```
unzip ParseBiosciences-Pipeline.1.3.1.zip; rm ParseBiosciences-Pipeline.1.3.1.zip
cd ParseBiosciences-Pipeline.1.3.1/
```
#### iv. Run install script in pipeline directory
```
bash ./install_dependencies_conda.sh -i -y
pip install --no-cache-dir ./
```
#### v. Check installation of Parse split-pipe. 
Calling the empty function should return the version.
```
split-pipe
```

## 5. Build reference genome
#### i. (Optional) Add GFP and mCherry to the human genome reference 
The ClonMapper Duo system is designed such that indexed sequences are included in the barcode information for each cell, for high confidence fusion cell assignment. While GFP and mCherry transcripts may not be picked up in every cell's transcriptome, GFP and mCherry sequences can be included in the annotation file so that these sequences can be mapped when detected. Mapping of these sequences will provide additional confidence in assignments of fusion cells. 

a. In the genomes folder, create a file called 'genome.gfp_mcherry.fa' containing sequence information for mCherry and eGFP
```
>MCHERRY
GAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATG
>EGFP
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA
```

b. In the genomes folder, create a file called 'genes.gfp_mcherry.gtf' containing annotation information for mCherry and eGFP (mCherry is 678 bp in length, GFP is 720 bp in length)
```
MCHERRY	unknown	gene	1	678	.	+	.	gene_id "MCHERRY"; gene_name "MCHERRY"; gene_biotype "protein_coding";
MCHERRY	unknown	transcript	1	678	.	+	.	gene_id "MCHERRY"; transcript_id "MCHERRY"; gene_name "MCHERRY"; gene_biotype "protein_coding";
MCHERRY	unknown	exon	1	678	.	+	.	gene_id "MCHERRY"; transcript_id "MCHERRY"; gene_name "MCHERRY"; gene_biotype "protein_coding";
EGFP	unknown	gene	1	720	.	+	.	gene_id "EGFP"; gene_name "EGFP"; gene_biotype "protein_coding";
EGFP	unknown	transcript	1	720	.	+	.	gene_id "EGFP"; transcript_id "EGFP"; gene_name "EGFP"; gene_biotype "protein_coding";
EGFP	unknown	exon	1	720	.	+	.	gene_id "EGFP"; transcript_id "EGFP"; gene_name "EGFP"; gene_biotype "protein_coding";
```

c. Concatenate the mCherry/GFP sequence and annotation files to the GRCh38 sequence and annotation files created new files called 'GRCh38_GFP_mCherry.genome.fa' and 'GRCh38_GFP_mCherry.genes.gtf'
```
gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz 
cat Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa genome.gfp_mcherry.fa > GRCh38_GFP_mCherry.genome.fa
gzip GRCh38_GFP_mCherry.genome.fa

gunzip Homo_sapiens.GRCh38.112.gtf
cat Homo_sapiens.GRCh38.112.gtf genes.gfp_mcherry.gtf > GRCh38_GFP_mCherry.genes.gtf
gzip GRCh38_GFP_mCherry.genes.gtf
```

#### ii. Build reference genome
Here, you will provide sequence and annotation files and a Parse compatible reference genome will be built. This step can take 30 minutes or more, launch screen or another program to ensure connection to server is not interrupted during execution.

a. If you added GFP and mCherry to the references, run this command from the 'parsealign' directory. Adding the flag '--dryrun' after the command will perform call checks without running.
```
split-pipe --mode mkref --genome_name hg38 --fasta ./genomes/GRCh38_GFP_mCherry.genome.fa.gz --genes ./genomes/GRCh38_GFP_mCherry.genes.gtf.gz --output_dir ./genomes/hg38_gfp_mch --nthreads 16
```

a. If you did not add GFP and mCherry to the reference, run this command to build the standard Hg38 human reference. Adding the flag '--dryrun' after the command will perform call checks without running.
```
split-pipe --mode mkref --genome_name hg38 --fasta ./genomes/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz --genes ./genomes/Homo_sapiens.GRCh38.112.gtf.gz --output_dir ./genomes/hg38
```

## 6. Run alignment for Parse Whole Transcriptome data

Most of this code and workflow is standard and adapted from Parse.  Each sublibrary alignment step can take 7 hours or more depending on your system. Launch screen or another program to ensure connection to server is not interrupted during execution, each sublibrary may take 7+ hours to complete. 
#### i. Make directory to store aligned results
```
mkdir aligned
```

#### ii. Create a space-delimited file called 'sample-info.txt' showing where each sample was loaded into the Parse stage 1 barcoding plate
```
H1806_initial A1-A4
H1806_1_presort A5-A7
H1806_2_presort A8-A10
H1806_3_presort A11-B1
H1806_4_presort B2-B4
H1806_1_ctrl B5
H1806_2_ctrl B6 
H1806_3_ctrl B7
H1806_4_ctrl B8
H1806_1_fusion B9
H1806_2_fusion B10
H1806_3_fusion B11
H1806_4_fusion B12
MB231_initial C1-C4
MB231_1_presort C5-C7
MB231_2_presort C8-C10
MB231_3_presort C11-D1
MB231_4_presort D2-D4
MB231_1_ctrl D5
MB231_2_ctrl D6
MB231_3_ctrl D7
MB231_4_ctrl D8
MB231_1_fusion D9 
MB231_2_fusion D10
MB231_3_fusion D11 
MB231_4_fusion D12
```

#### iii. Perform alignment.
a. Each sublibrary can be run one at a time. Adding the flag '--dryrun' will perform checks without running. This example shows the call for a single Parse sublibrary, "WT-1".
```
conda activate spipe

split-pipe \
--mode all \
--chemistry v2 \
--genome_dir ./genomes/hg38_gfp_mch/ \
--fq1 ./expdata/WT/WT-1.R1.raw.fastq.gz \
--fq2 ./expdata/WT/WT-1.R2.raw.fastq.gz \
--output_dir ./aligned/WT-1 \
--nthreads 16 \
--samp_list sample-list.txt
```

b. Or a script can be used to run each sample in series:
```
#!/bin/bash
# this script performs alignment for each Parse sublibrary one at a time, storing the output of each 

samples=("1" "2" "3" "4" "6" "7" "8")

for i in "${samples[@]}"; do
	# find files for each Parse sublibrary and set output dir to sample name
	fq1="./expdata/WT/WT-${i}.R1.raw.fastq.gz"
	fq2="./expdata/WT/WT-${i}.R2.raw.fastq.gz"
	output_dir="./aligned/WT-${i}"

	# run Parse alignment
	split-pipe \
	--mode all \
	--chemistry v2 \
	--genome_dir ./genomes/hg38_gfp_mch/ \
	--fq1 "$fq1" \
	--fq2 "$fq2" \
	--output_dir "$output_dir" \
	--nthreads 8 \
	--samp_list sample-list.txt
done
```
#### iii. Combine results from each sublibrary
```
split-pipe --mode comb \
--sublibraries \
./aligned/WT-1 \
./aligned/WT-2 \
./aligned/WT-3 \
./aligned/WT-4 \
./aligned/WT-6 \
./aligned/WT-7 \
./aligned/WT-8 \
--output_dir ./aligned/combined_WT \
--nthreads 16
```


## 7. Run alignment for Parse CRISPR Detect data

The CRISPR detect pipeline requires a guide RNA list to be provided. We obtain a list of barcodes from targeted barcode sequencing (See Step 2) and prepped a csv in the following form (see step 3)

`formatted_guide_list.csv`
```
#Guide_Name,Prefix,Guide_Sequence,Suffix,Target_Gene
ATCGCCAGGGCGAGGGGTAG,ATCTTGTGGAAAGGACGAAACACCG,ATCGCCAGGGCGAGGGGTAG,GTTTTAGAGCTAGAAATAGCAAGTT,ATCGCCAGGGCGAGGGGTAG
GCTGTGGTAACAGGGCACGG,ATCTTGTGGAAAGGACGAAACACCG,GCTGTGGTAACAGGGCACGG,GTTTTAGAGCTAGAAATAGCAAGTT,GCTGTGGTAACAGGGCACGG
```

The 'Prefix' (the 5' front flanking region of the barcode for all barcodes) is: ATCTTGTGGAAAGGACGAAACACCG
The suffix (the 3' rear flanking region of the barcode for all barcodes) is: GTTTTAGAGCTAGAAATAGCAAGTT

The CRISPR detect pipeline requires a guide RNA list to be provided. We obtain a list of barcodes from targeted barcode sequencing (See Step 2) and need to prep a csv in the following form: 

`formatted_guide_list.csv`
```
#Guide_Name,Prefix,Guide_Sequence,Suffix,Target_Gene
ATCGCCAGGGCGAGGGGTAG,ATCTTGTGGAAAGGACGAAACACCG,ATCGCCAGGGCGAGGGGTAG,GTTTTAGAGCTAGAAATAGCAAGTT,ATCGCCAGGGCGAGGGGTAG
GCTGTGGTAACAGGGCACGG,ATCTTGTGGAAAGGACGAAACACCG,GCTGTGGTAACAGGGCACGG,GTTTTAGAGCTAGAAATAGCAAGTT,GCTGTGGTAACAGGGCACGG
```

Here, we set the barcode sequence to be the 'Guide_Name', 'Guide_Sequence', and 'Target_Gene'
The 'Prefix' (the 5' front flanking region of the barcode for all barcodes) is: ATCTTGTGGAAAGGACGAAACACCG
The suffix (the 3' rear flanking region of the barcode for all barcodes) is: GTTTTAGAGCTAGAAATAGCAAGTT

#### i. Add # to the header line to our formatted_guide_list.csv to match Parse formatting
```
sed '1s/^/#/' formatted_guide_list.csv > formatted_guide_list_hashed.csv
```

Our 'formatted_guide_list_hashed.csv' file should now look like this:
```
#Guide_Name,Prefix,Guide_Sequence,Suffix,Target_Gene
g001_h1806,ATCTTGTGGAAAGGACGAAACACCG,ATCGCGTATCCAGAAGTACG,GTTTTAGAGCTAGAAATAGCAAGTT,ATCGCGTATCCAGAAGTACG
g002_h1806,ATCTTGTGGAAAGGACGAAACACCG,GACAACGCTGCGGGGAGGGG,GTTTTAGAGCTAGAAATAGCAAGTT,GACAACGCTGCGGGGAGGGG
g003_h1806,ATCTTGTGGAAAGGACGAAACACCG,GACAACAATTTTGTCACATT,GTTTTAGAGCTAGAAATAGCAAGTT,GACAACAATTTTGTCACATT
g004_h1806,ATCTTGTGGAAAGGACGAAACACCG,GCTGTCCGTCTTTGTGGTAG,GTTTTAGAGCTAGAAATAGCAAGTT,GCTGTCCGTCTTTGTGGTAG
g005_h1806,ATCTTGTGGAAAGGACGAAACACCG,ACGCAGGCACGACGAGCCTA,GTTTTAGAGCTAGAAATAGCAAGTT,ACGCAGGCACGACGAGCCTA
```

#### ii. Move back to your main directory ('parsealign') and perform a dry run to ensure there are no errors in the CRISPR Detect pipeline inputs
```
conda activate spipe

split-pipe --mode all \
--chemistry v2 \
--crispr \
--crsp_guides ./pycash_analysis_clean/formatted_guide_list_hashed.csv \
--parent_dir ./aligned/WT-1 \
--output_dir ./aligned/CR-1 \
--fq1 ./expdata/CRISPRdetect/CR-1.R1.raw.fastq.gz \
--fq2 ./expdata/CRISPRdetect/CR-1.R2.raw.fastq.gz \
--nthreads 16 \
--crsp_tscp_thresh 3 \
--dryrun
```

#### iii. Run CRISPR Detect alignment 
If the dry run was successful, launch screen and a script to run all samples 

`run_crispr_align.sh`
```
#!/bin/bash
# this script performs alignment for each Parse CRISPR Detect sublibrary one at a time, storing the output of each 

samples=("1" "2" "3" "4" "6" "7" "8")

for i in "${samples[@]}"; do
	# run Parse alignment
	split-pipe --mode all \
	--chemistry v2 \
	--crispr \
	--crsp_guides ./pycash_analysis_clean/formatted_guide_list_hashed.csv \
	--parent_dir ./aligned/WT-${i} \
	--output_dir ./aligned/CR-${i} \
	--fq1 ./expdata/CRISPRdetect/CR-${i}.R1.raw.fastq.gz \
	--fq2 ./expdata/CRISPRdetect/CR-${i}.R2.raw.fastq.gz \
	--crsp_read_thresh 3 \
	--crsp_tscp_thresh 3 \
	--nthreads 16
done
```

```
bash ./scripts/run_crispr_align.sh
```

#### vi. Combine CRISPR detect outputs 
```
split-pipe --mode comb \
--parent_dir ./aligned/combined_WT \
--output_dir ./aligned/combined_CR \
--sublibraries \
./aligned/CR-1 \
./aligned/CR-2 \
./aligned/CR-3 \
./aligned/CR-4 \
./aligned/CR-6 \
./aligned/CR-7 \
./aligned/CR-8 \
--nthreads 12
```
