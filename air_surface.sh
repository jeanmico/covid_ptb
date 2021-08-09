#! /usr/bin/env bash
#SBATCH --job-name=air-surface
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out

/software/c4/wittelab/containers/geospatial-3.6.2/rocker_r-geospatial.img Rscript air_surface.R

[[ -n "$SLURM_JOB_ID" ]] && sstat --format="JobID,AveCPU,MaxRSS,MaxPages,MaxDiskRead,MaxDiskWrite" -j "$SLURM_JOB_ID"
