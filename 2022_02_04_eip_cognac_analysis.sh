#!/bin/sh
#SBATCH --job-name=2022_02_04_eip_cognac_analysis
#SBATCH --mail-user=rcrawfo@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=16g
#SBATCH --time=12-00:00:00
#SBATCH --account=esnitkin1
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/umms-esnitkin/ryan/eip_hgt_analysis/lib/%x-%j.log

cd /nfs/turbo/umms-esnitkin/ryan/eip_hgt_analysis/lib/
Rscript 2022_02_04_eip_cognac_analysis.R
