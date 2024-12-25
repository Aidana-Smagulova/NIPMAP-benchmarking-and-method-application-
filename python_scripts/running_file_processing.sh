#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mem=50gb
#SBATCH --export=NONE

SCRIPT_PATH="/gpfs/bwfor/work/ws/hd_vb248-nipmap/NIPMAP/processing_files.py"
python3 $SCRIPT_PATH
