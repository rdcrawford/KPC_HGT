#!/bin/sh
#SBATCH --job-name=2022_02_04_run_fastTree_cognac_nt
#SBATCH --mail-user=rcrawfo@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=750gb
#SBATCH --time=12-00:00:00
#SBATCH --account=esnitkin1
#SBATCH --partition=largemem
#SBATCH --output=/nfs/turbo/umms-esnitkin/ryan/eip_hgt_analysis/lib/%x-%j.log

cd /nfs/turbo/umms-esnitkin/ryan/eip_hgt_analysis/analysis/2022_02_04_eip_cognac_analysis/
fasttree -nt concatenated_gene_nt_alignment.fasta > concatenated_gene_nt_alignment_fasttree.tre
