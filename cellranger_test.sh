#bin/bash
#BSUB -J cellranger_test2
#BSUB -q normal
#BSUB -n 2
#BSUB -R "rusage[mem=4000]"
#BSUB -W 24:00
#BSUB -N
#BSUB -o /home/r22101006/demuxly/cellranger_test2.out
#BSUB -e /home/r22101006/demuxly/cellranger_test2.err

cd /home/r22101006/demuxly

cellranger count \
	--id=cellranger_test2 \
	--create-bam=true \
	--transcriptome=/home/r22101006/demuxly/reference/refdata-gex-GRCh38-2024-A \
	--fastqs=/home/r22101006/demuxly/data/pbmc_1k_v3_fastqs \
	--sample=pbmc_1k_v3
