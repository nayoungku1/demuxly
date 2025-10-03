# DEMUXLY 🧶: Single cell RNA-seq data demultiplexing algorithm



## Project Structure
```
demuxly/
├── data/
│   ├── pbmc_1k_v3_fastq/
│   │   ├── pbmc_1k_v3_S1_L001_I1_001.fastq.gz
│   │   ├── pbmc_1k_v3_S1_L001_R1_001.fastq.gz
│   │   ├── ...
│   │   └── pbmc_1k_v3_S1_L002_R2_001.fastq.gz
│   └── sample_data/
│       └── small_test.fa
├── reference/
│   ├── refdata-gex-GRCh38-2024-A/ - human reference genomes
│   │   ├── fasta/
│   │   │   ├── genome.fa
│   │   │   └── genome.fai
│   │   ├── genes/
|   │   │   └──  genes.gtf
│   │   └── star/
│   └── STAR_index/
├── results/
│   └── output/
├── yard/
│   └── apps/
│       ├── cellranger-6.0.2/
│       ├── cellranger-7.0.0/
│       ├── cellranger-8.0.1/
│       ├── cellranger-9.0.1/
│       └── check_install/
└── README.md
```
## Installation / Setup
* Anaconda/Miniconda enviornment:
```bash
conda env create -f environment.yml
conda activate demuxly
```
* cellranger:
    ```bash
    mkdir yard/apps
    cd yard/apps
    ```
  * `version 9.0.1`
    ```bash
    wget -O cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1759546977&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=UC2wRb8nMupMap9f9cWBqD~BQCtUPhF6TVuWchIMA9lwlG71eDnYTKEKUa6BpbO6TiQRqqv3YbQ3V7whx7wrw4ws~oPmsXcc~ugB3H3Lajso2Lf8QM99KZL6gU25GMd4xQs~F33bALSp6LaEuLxbBsIwTb54SpMNrtbcjrkM3cT2jMCiGxIbJ7hEtd4cRA3rnTC~MwwsiNMta5QeeDGJnjtYT~nXzQPf6UkcMAn4TCvPlcl70NxI87tEiEm~SvIqAjO4b3K8GeP4H3RQrDyU6o-XdajYK5sWQPwSWHbwkZGqOQ9lZ4j~8ofNstIm03k2gqwvMHkKvSuKhgAsCBIA7g__"
    tar -xzvf cellranger-9.0.1.tar.gz

    ```
  * `version 8.0.1`
    ```bash
    curl -o cellranger-8.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.1.tar.gz?Expires=1759546978&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=BU6-nro5ZJbBk7ARtOpnLWPeEDzDxZDeze-KKLbxpcPA5cDB9QHCMODPdFF2MoxPveoLGzxQfdzQarMw3P1zIlghHgzkng4SzFWd9NDqJc-jEkSuNjDG3dkZjmebCqODJVDwufW3W9RJGj8KTH-5PweSi93x0uiSFW67eUgJDF72ZTpFy3fJRMtgHb4ycqrmOBvPevwWNWTjyg1k8aR-DAAmM9lLFeI9Dz21ITjQRexp8yXXoWQ-Witye0x1IvvmiL0K8nDT78CPdhuFghF1sC~XiDPIDCAMW23s7~jTwyHX4BP-6RdEHTOLG5WGxB04XtT4xr9B4HK0tvlIyuWupg__"
    tar -xzvf cellranger-8.0.1.tar.gz
    ```
  * `version 7.0.0`
    ```bash
    wget -O cellranger-7.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.0.0.tar.gz?Expires=1759546978&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=ZSQ5sqOUHi98htHiQaLpVMUdT3VINGTmC4kkjPeRt2F0od4cgb55H7QdOZEV9OFGVLpQrHLqnlUxexWZveWTxSI-8n15t44KMPlQNM0PsMRdRqduar1BYhv0jvZzumMNxMeEa6JlsrK8yxWP7qfCWYgvDG7MGf9UdV697Bov-Ci6rh0EQewEOV4GFEwWCzMPwqK0KgtD7owFFVXsiEQUa6mtNMt2vbOEgFVtaRT5TbNGTcSKjLSRsI~zbIt5iRSvcpNrpjuJ7snpibhdxh3BHe7~2dBAgEki4LDMXzr1pjEODZy369uDDvSCJlwoEpODRawFOs97hN6Xa6fFTa0N5g__"
    tar -xzvf cellranger-7.0.0.tar.gz
    ``` 
  * `version 6.0.2`
    ```bash
    wget -O cellranger-6.0.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.0.2.tar.gz?Expires=1759546978&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=U1v4UX6ur0EYEDTz-ML98J2GBW4GwfAOxpN8hUeWTNqAScDouhi2CvzapVPP~pCa6vT9NTy6d3DRDGkXD4aHhF4cFUTIuu3UAoBLNuQuwRRm4eLOMFfmLhKK7rW5N9pWERCcj2UOx2LuJR8y7GtAFPDyhlnmXYXsx3tog0DJP7K2dvVqDj~hIOr~8SVkbnX7okAvAjVQ7mNBqFY1kxv5a-e65tKc0fxlGoW0lu-t7fvf2xHfcPRm5QCjEFi-gJA479NRfJM8VWyl4Hh6THtZCfvbOpcS79mU0T~~C-gYRsOTBtuZa-fZmaQd8GdTLGQsGQC4Ym1ioeYmjLVW5R0f7w__"
    tar -xzvf cellranger-6.0.2.tar.gz
    ```
## Sample Dataset
* **Human Reference Genome** (GRCh38) - 2024-A:  
```bash
wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz
gunzip refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz
```
* **STAR genome index** for Alignment using STARsolo
```bash
STAR --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir STAR_index \
  --genomeFastaFiles refdata-gex-GRCh38-2020-A/fasta/genome.fa \
  --sjdbGTFfile refdata-gex-GRCh38-2020-A/genes/genes.gtf \
  --sjdbOverhang 90
```  
* **1k PBMCs from a Healthy Donor (v3 chemistry)**: Universal 3' Gene Expression dataset analyzed using Cell Ranger 3.0.0; single sample
  * 1,222 cells detected
  * Sequenced on Illumina NovaSeq with approximately 54,000 reads per cell
  * 28bp read1 (16bp Chromium barcode and 12bp UMI), 91bp read2 (transcript), and 8bp I7 sample barcode
  * run with --expect-cells=1000 
  * You can manually download the file from [here](https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0)
```bash
# FASTQ
mkdir data
cd data

wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
tar -xvf pbmc_1k_v3_fastqs.tar

# Feature / Cell Matrix (.h5, raw) using Cell Ranger
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_raw_feature_bc_matrix.h5
```

* **10k Human PBMCs Multiplexed, 2 CMOs**: Universal 3' Gene Expression dataset; **MULTIPLEXED DATA**
  * run with the expect number of cells as 10000 
  * You can manually download the file from [here](https://www.10xgenomics.com/datasets/10-k-human-pbm-cs-multiplexed-2-cm-os-3-1-standard-6-0-0)
```bash
# FASTQ

# Input Files
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex_fastqs.tar
wget https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex_config.csv
wget https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex_count_feature_reference.csv

tar -xvf SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex_fastqs.tar
```
  * To run with ```cellranger-8.0.1``` or above, you should add a value for `create-bam` in the configuration file's `[gene-expression]`:
```
vim SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex_config.csv
...
[gene-expression]
reference,/path/to/demuxly/reference/refdata-gex-GRCh38-2024-A
expect-cells,10000
create-bam,true
...
```      

* **10k PBMCs from a Healthy Donor (v3 chemistry)**: Universal 3' Gene Expression dataset analyzed using Cell Ranger 3.0.0
  * Peripheral blood mononuclear cells (PBMCs) from a healthy donor.
  * PBMCs are primary cells with relatively small amounts of RNA (~1pg RNA/cell). 
  * You can manually download this BAM file from [here](https://www.10xgenomics.com/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0)  
```bash
  
  cd data
  wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam
  ```
  
