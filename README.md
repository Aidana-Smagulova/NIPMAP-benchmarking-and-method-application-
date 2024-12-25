# NIPMAP

Lukas data: 
regionprops_annotated - csv files with cell coordinates in images, 1000x1000 px, 1000x1000µm
intensities - csv files with marker intensities in images 

scripts to process data: nipnap_working.ipynb or interactive.sh with regionprops_csv.py 

script to run 4 niches with 1/4 of patients: 4n_nipmap.py 

py_wrapper_nipmap.r - R script to process data, more updated version is stored locally 

process: 
- take updated csv's - two folders (intensities & regionsprops) with csv files for each patient, sample and ROI. names in format: TS-373_IMC78_UB_001.csv (TS-373_IMC***_***_00*.csv)


parameters to change: 
- protein markers
- cell types (phenotypes)
- number of niches
- radii of sites
- number of sites in the script 
