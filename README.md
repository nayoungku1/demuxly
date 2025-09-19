# DEMUXLY 🧶: Single cell RNA-seq data demultiplexing algorithm

## Setup
```bash
conda env create -f environment.yml
```

## Sample Dataset
* **10k PBMCs from a Healthy Donor (v3 chemistry)**: Universal 3' Gene Expression dataset analyzed using Cell Ranger 3.0.0
  * Peripheral blood mononuclear cells (PBMCs) from a healthy donor.
  * PBMCs are primary cells with relatively small amounts of RNA (~1pg RNA/cell). 
  ```bash
  wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_possorted_genome_bam.bam
  ```
  * You can manually download this BAM file from [here](https://www.10xgenomics.com/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0)
