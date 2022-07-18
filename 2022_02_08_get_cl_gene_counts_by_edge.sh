#!/bin/sh
#SBATCH --job-name=2022_02_08_get_cl_gene_counts_by_edge
#SBATCH --mail-user=rcrawfo@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=24gb
#SBATCH --time=12-00:00:00
#SBATCH --account=esnitkin1
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/umms-esnitkin/ryan/eip_hgt_analysis/lib/%x-%j.log

cd /nfs/turbo/umms-esnitkin/ryan/eip_hgt_analysis/lib/
Rscript 2022_02_08_get_cl_gene_counts_by_edge.R
#sbatch 2022_02_09_find_KPC_acq_enriched_edges.sh
#sbatch 2022_02_11_find_cluster_precision.sh
