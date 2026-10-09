#!/usr/bin/bash -l
#SBATCH --job-name=snake_figs_sub
#SBATCH --output=snake_log/job_%j_%x.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90
#SBATCH --mail-user=mter@bio.aau.dk

# Exit on first error, unset variable, or failed pipe
set -euo pipefail

mkdir -p snake_log

# Any env with Snakemake >= 9 and the SLURM executor plugin works (the one
# made from Snakemake_env.yml is snakemake_figs); change the name to yours.
conda activate snakemake_figs

# Snakemake itself needs almost nothing: it submits each rule as its own
# SLURM job via the executor plugin. 1 CPU + 4G is plenty.

# No setup step needed -- each rule downloads its own container the first
# time it runs (see README "Containers").
snakemake --profile profile

# If the workflow crashed and the directory is locked, run once with --unlock,
# then remove the flag and rerun:
#snakemake --profile profile --unlock
