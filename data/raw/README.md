# Raw Data

The raw data used in this project were obtained from The Cancer Genome Atlas (TCGA) through the Genomic Data Commons (GDC) Data Portal.

The following files are required to reproduce the analyses:

- TCGA-BRCA.star_tpm.tsv.gz
- TCGA-BRCA.clinical.tsv.gz
- TCGA-BRCA.survival.tsv.gz
- gencode.v36.annotation.gtf.gene.probemap

These files are not included in this repository because they are large and publicly available.

## Download

The TCGA BRCA expression, clinical, and survival data can be downloaded from:

https://portal.gdc.cancer.gov/

The gene annotation file can be downloaded from:

https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap

After downloading, place all files inside this directory:

data/raw/