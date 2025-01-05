# NIPMAP

The data: 
regionprops_annotated - csv files with cell coordinates in images, 1000x1000 px, 1000x1000µm
intensities - csv files with marker intensities in images 

script to run 4 niches with 1/4 of patients: 4n_nipmap.py 

py_wrapper_nipmap.r - R script to process data, more updated version is stored locally 

process: 
- produce csvs from the ann data object: two folders (intensities & regionsprops) with csv files for each patient, sample and ROI. names in format: TS-373_IMC78_UB_001.csv (TS-373_IMC*patient ID*_*disease*_*ROI*.csv)
- intensities folder contains files for each patient ID with phenotypic markers and marker intensities for each cell (object).
- regionsprops folder contains files for each patient, with X and Y coordinates and phenotypes for each cell (object)
- the script preprocessing.py will rename columns as in the NIPMAP tutorial (https://github.com/jhausserlab/NIPMAP/blob/main/README.md#inputs) and create a main 

parameters to change: 
- protein markers
- cell types (phenotypes)
- number of niches
- radii of sites
- number of sites in the script 
