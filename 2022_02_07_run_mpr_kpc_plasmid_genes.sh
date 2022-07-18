#!/bin/sh
#SBATCH --job-name=2022_02_07_run_mpr_kpc_plasmid_genes
#SBATCH --mail-user=rcrawfo@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=6gb
#SBATCH --time=1-00:00:00
#SBATCH --account=esnitkin1
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/umms-esnitkin/ryan/eip_hgt_analysis/lib/%x-%j.log

cd /nfs/turbo/umms-esnitkin/ryan/eip_hgt_analysis/lib
Rscript 2022_02_07_run_mpr_kpc_plasmid_genes.R
#Rscript 2022_02_07_find_best_hclust_height_threshold.R
#sbatch 2022_02_08_get_cl_gene_counts_by_edge.sh
