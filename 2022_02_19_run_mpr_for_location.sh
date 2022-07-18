#!/bin/sh
#SBATCH --job-name=2022_02_19_run_mpr_for_location
#SBATCH --mail-user=rcrawfo@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=6gb
#SBATCH --time=00-01:00:00
#SBATCH --account=esnitkin1
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/umms-esnitkin/ryan/eip_hgt_analysis/lib/%x-%j.log

cd /nfs/turbo/umms-esnitkin/ryan/eip_hgt_analysis/lib
Rscript 2022_02_19_run_mpr_for_location.R
