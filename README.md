# DEMUXLY 🧶: Single cell RNA-seq data demultiplexing algorithm

## Setup
```bash
conda env create -f environment.yml
conda activate demuxly
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
* **1k PBMCs from a Healthy Donor (v3 chemistry)**: Universal 3' Gene Expression dataset analyzed using Cell Ranger 3.0.0
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
* **10k PBMCs from a Healthy Donor (v3 chemistry)**: Universal 3' Gene Expression dataset analyzed using Cell Ranger 3.0.0
  * Peripheral blood mononuclear cells (PBMCs) from a healthy donor.
  * PBMCs are primary cells with relatively small amounts of RNA (~1pg RNA/cell). 
  * You can manually download this BAM file from [here](https://www.10xgenomics.com/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0)  
```bash
  
  cd data
  wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam
  ```
  
