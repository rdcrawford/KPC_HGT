#!/bin/sh
#SBATCH --job-name=2022_02_06_get_transition_edges_from_cognac_tree_mpr
#SBATCH --mail-user=rcrawfo@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=24gb
#SBATCH --time=12-00:00:00
#SBATCH --account=esnitkin1
#SBATCH --partition=largemem
#SBATCH --output=/nfs/turbo/umms-esnitkin/ryan/eip_hgt_analysis/lib/%x-%j.log

cd /nfs/turbo/umms-esnitkin/ryan/eip_hgt_analysis/lib
Rscript 2022_02_06_get_transition_edges_from_cognac_tree_mpr.R
