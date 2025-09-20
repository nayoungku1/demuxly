#!/bin/bash
#BSUB -J barcode_rank         
#BSUB -q normal                
#BSUB -n 1                    
#BSUB -R "rusage[mem=2000]"    
#BSUB -W 7:00                 
#BSUB -o /home/r22101006/demuxly/barcode_rank/barcode_rank.out        
#BSUB -e /home/r22101006/demuxly/barcode_rank/barcode_rank.err         
#BSUB -N

cd /home/r22101006/demuxly

conda run -n demuxly python barcode_ranking.py \
    --fastq_folder data/pbmc_1k_v3_fastqs \
    --cellranger_h5 data/pbmc_1k_v3_raw_feature_bc_matrix.h5 \
    --output_dir barcode_rank \
    --manual_topk 1000 \
    --expect_cells 1000
