# README

This document aims to summarise the result of a 2-month computational internship. The contents are:
- NIPMAP results reproduction with MIBI images of TNBC of Keren et al (2021)
- re-application of NIPMAP on the bone cancer dataset (from Lukas Hatscher) 
- comparison and interpretation of the method output & conclusions

This README file bases on the original [NIPMAP GitHub repository](https://github.com/jhausserlab/NIPMAP/tree/main) and [publication](https://www.nature.com/articles/s41467-023-42878-z) by El Marrahi et al. 2023


## Method description 
NIPMAP - NIche Phenotype MAPping (NIPMAP) analysis from spatial multiplex data.

Tissue microenvironments play a significant role in pathology. 

NIPMAP is a spatial mathod that works to discover niches (clusters of cells with similar profile) using principles from the community ecology approach, PCA and archetypes analysis. 
For this method, each cell type is treated like a species in an oncological niche, where its abundance can be inferred. Based on thus abundance, analysis like PCA can be performed. The PCA space is then fitted onto a simpled figure (using AA), where simplex tips are archetypes (extreme cases of niches). Each cell will find itself somewhere on the simplex. Cells right at the tips will most probably belong to that archetype and be situated in the corresponding niche. However, cells that would lie on the edges of the simplex, more specifically between two archetypes, most likely lay between two niches and will be, therefore, assigned to an interface between the two niches. 

## index? 

## Requirements 

Tested on: OS: Red Hat Enterprise Linux 8.10 (Ootpa)

Use the minimal_env.yml file to re-create working conda environment for NIPMAP. 

Or run: 

```bash
    conda create --name NIPMAP python=3.7.13
```
and then:

```bash
    pip install matplotlib
    pip install scipy
    pip install pandas==1.3.5
    pip install numpy>=1.17.3
    pip install scikit-learn
    pip install seaborn
    pip install qpsolvers==1.9.0
``` 


Further description of all packages in the environment can be found in the requirements.txt file. 

* R libraries:
    ```
    - R 4.1.3
    - RStudio
    run: 
    pkgs <- c("tidyverse","ggplot2","ade4","factoextra","plotly","igraph","reshape2","ggrepel","viridis","fdrtool","pheatmap","cluster","broom","pROC","ggpubr","devtools","ggridges")
    install.packages(pkgs)
    ```

## Proof of concept: 
ran NIPMAP on the original data with the GITHUB repo cloned: 
maybe compare the results to the original ones? i dunno. have to run it again i guess, cause my dumbass forgot to save any outputs... 

Data: Multiplex Ion Imaging on 41 Triple Negative Breast tumors from [Keren et al, Cell(2018)](10.1016/j.cell.2018.08.039) and In Situ Sequencing data on human lung development from [Sountoulidis et al, 2022](https://doi.org/10.1101/2022.01.11.475631)

## Re-application of the method on Lukas' data 

### Data: 
regionprops_annotated - csv files with cell coordinates in images, 1000x1000 px, 1000x1000µm
intensities - csv files with marker intensities in images 

```
.
├── 7n_output
├── TMENS_analysis
│   └── data
|       ├── cell_data.csv #big csv with all of the cells & marker intensities 
│       └── cell_positions_data #contains all patients' files with data on X and Y coordinates of each cell + the renaming csv 
└── updated_csvs #raw data 
    ├── intensities #cell phenotypic markers & intensities 
    └── regionprops #cell coordinates, phenotype 
```


process: 
- produce csvs from the ann data object: two folders (intensities & regionsprops) with csv files for each patient, sample and ROI. names in format: TS-373_IMC78_UB_001.csv (TS-373_IMC*patient ID*_*disease*_*ROI*.csv)
- intensities folder contains files for each patient ID with phenotypic markers and marker intensities for each cell (object).
- regionsprops folder contains files for each patient, with X and Y coordinates and phenotypes for each cell (object)
- the script preprocessing.py will rename columns as in the NIPMAP tutorial (https://github.com/jhausserlab/NIPMAP/blob/main/README.md#inputs) and create a main csv where names of the files are coded into integers 

1. One CSV file per sample named *patient\<Patient ID\>_cell_positions.csv*, with cells as rows and the following columns:
* *x*, x position of the cell in the image
* *y*, y position of the cell in the image
* *label*, an integer index for the cell, to cross reference with the *cellData.csv* file (below)
* *cell_type*, a string describing the cell type
2. One CSV file named *cellData.csv* with cells as rows and the following columns:
* *cellLabelInImage*, an integer index for the cell, to cross reference with the *patient\<Patient ID\>_cell_positions.csv* files (above)
* *SampleID*, the sample of origin, matching the *\<Patient ID\>* of the *patient\<Patient ID\>_cell_positions.csv* files
* one column per marker containing the intensity of each marker (for example MHCI, CD68, CD45RO, ...)
The repository contains example CSV files from the Keren et al. and Sountoulidis et al. studies.


### Running the analysis 

bash running_nipmap.sh - submitting a batch job 
uses the python script - running_nipmap.py

### parameters to choose: 
parameters to change: 
- protein markers
- cell types (phenotypes)
- number of niches
- radii of sites
- number of sites in the script

### Working with nipmap 
The python script will generate these outputs as json files:
* Generation of sites and cellular abundances within them (+ radius size selection)
* PCA and Archetype Analysis
* Niche identification
* Niche segmentation of images


### outputs 
NOTEBOOKS 

2. Open *py_wrapper_nipmap.r* script. Set the parameters in the script header. Execute the script in Rstudio or on the command line:
  ```bash
  R --vanilla < py_wrapper_nipmap.r
  ```
  This script produces these outputs including figures: 
* Niche-phenotype mappping

i just took one json file to turn it into a giant csv and work further with it 


Note: one sample = one image
Note #2: NIPMAP doesn't aim to correct cell segmentation error or cell type mis-assignments, this needs to be addressed prior to niche-phenotype mapping. 


### results / interpretation / scatterplots 
COMPARING WITH K MEANS CLUSTERING 
THE INTERFACE ISSUE 
ASSIGNING NICHES TO CELLS 


El Marrahi, A., Lipreri, F., Kang, Z. et al. NIPMAP: niche-phenotype mapping of multiplex histology data by community ecology. Nat Commun 14, 7182 (2023). https://doi.org/10.1038/s41467-023-42878-z
