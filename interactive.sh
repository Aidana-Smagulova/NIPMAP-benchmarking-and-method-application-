#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --mem=50gb
#SBATCH --export=NONE

SCRIPT_PATH="/gpfs/bwfor/work/ws/hd_vb248-nipmap/NIPMAP/scatterplot.py"
python3 $SCRIPT_PATH
