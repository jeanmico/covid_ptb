#! /usr/bin/env bash
#SBATCH --job-name=extract-values-any
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out

/software/c4/wittelab/containers/geospatial-3.6.2/rocker_r-geospatial.img Rscript extract_values_any.R

[[ -n "$SLURM_JOB_ID" ]] && sstat --format="JobID,AveCPU,MaxRSS,MaxPages,MaxDiskRead,MaxDiskWrite" -j "$SLURM_JOB_ID"
