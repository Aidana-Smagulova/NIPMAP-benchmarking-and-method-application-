#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem=1000gb
#SBATCH --export=NONE

source /home/hd/hd_hd/hd_vb248/.conda/envs/nipmap/activate
SCRIPT_PATH="/gpfs/bwfor/work/ws/hd_vb248-nipmap/NIPMAP/7niches_nipmap.py"
python3 $SCRIPT_PATH
deactivate
