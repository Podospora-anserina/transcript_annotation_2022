#!/bin/bash

#SBATCH --mem=50G
#SBATCH --cpus-per-task=10



set -euo pipefail


module load bedtools/2.30.0

work_dir="/shared/projects/minomics/pgrognet/RNAseq_Podo"

#script filter_sup_1500.R to get bed files with only nTARs bigger than 1.5kb

srun bedtools getfasta -fi /shared/projects/minomics/pgrognet/GenomePodo/genomePodoMatPlus.fasta -bed $work_dir/06-SCT_results/annotations/nTARs_gene_prediction/nTARs_1500.bed -fo $work_dir/06-SCT_results/annotations/nTARs_gene_prediction/nTARs_1500.fa


grep -v "^>" $work_dir/06-SCT_results/annotations/nTARs_gene_prediction/nTARs_1500.fa | awk 'BEGIN { ORS=""; print ">All_nTARs\n" } { print }' > $work_dir/06-SCT_results/annotations/nTARs_gene_prediction/nTARs_1500_merged.fasta