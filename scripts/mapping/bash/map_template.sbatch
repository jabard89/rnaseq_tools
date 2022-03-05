#!/bin/bash
#SBATCH --job-name=JB${BASHINDEX}
#SBATCH --output=/home/jbard/beagle3-dadrummond/jbard/211216/STAR/scripts/sbatch/map_JB${BASHINDEX}.out
#SBATCH --error=/home/jbard/beagle3-dadrummond/jbard/211216/STAR/scripts/sbatch/map_JB${BASHINDEX}.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --account=pi-dadrummond
#SBATCH --mem=30000

bash /home/jbard/beagle3-dadrummond/jbard/211216/STAR/scripts/sbatch/220303_star_JB${BASHINDEX}.sh ${BASHINDEX}
