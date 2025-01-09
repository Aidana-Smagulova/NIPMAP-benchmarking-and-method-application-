# README

This document aims to summarise the result of a 2-month computational internship. The contents include:
- application of NIPMAP on the bone cancer dataset (from Lukas Hatscher) 
- comparison and interpretation of the method output & conclusions

This README file bases on the original [NIPMAP GitHub repository](https://github.com/jhausserlab/NIPMAP/tree/main) and [publication](https://www.nature.com/articles/s41467-023-42878-z) by El Marrahi et al. 2023

## Index 


## Method description 
NIPMAP - NIche Phenotype MAPping (NIPMAP) analysis from spatial multiplex data.

Tissue microenvironments play a significant role in pathology. 

NIPMAP is a spatial mathod that works to discover niches (clusters of cells with similar profile) using principles from the community ecology approach, PCA and archetypes analysis. 
For this method, each cell type is treated like a species in an ecological niche. The whole image is sampled randomly and uniformly in small circular cites (radius can be decided manually), and abundancies of each cell types are calculated inside of the sampling cites. Based on this, an abundance matrix of cell types in sampling sites is created. 
This abundance matrix is then used to perform PCA (paper shows that 3 components are enough to capture the variance). After that, the PCA space is fitted onto a simplex figure using archetypal analysis (AA). An archetype in this case is a niche extreme case, when only one cell profile is abundant in a niche to 100% (say, only cancerous or only immune cells). Therefore, simplex tips are archetypes. If a sampling site falls on a samplex tip, the archetype weight will be 100% for this exact site. However, if a sampling site falls on an edge of the figure, it will not be allocated to a single archetype with 100% certainly but rather will represent a mixture of two near-lying archetype weights in percent. Sampling sites right in between of two archetypes (50% for each of them) are determined to be interfaces - regions that lay between any of the two niches. 

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

## Method application

### Data: 
```
.
├── 7n_output
├── TMENS_analysis
│   └── data
|       ├── cellData.csv #big csv with all of the cells & marker intensities 
│       └── cell_positions_data #contains all patients' files with data on X and Y coordinates of each cell + the renaming csv 
└── updated_csvs #raw data 
    ├── intensities #cell phenotypic markers & intensities 
    └── regionprops #cell coordinates, phenotype 
```

Data preprocessing steps: 
- produce csvs from the ann data object: two folders (intensities & regionsprops) with csv files for each patient, sample and ROI. names in format: TS-373_IMC78_UB_001.csv (TS-373_IMC_patientID_disease_ROI.csv)
- *intensities* folder contains files for each patient ID with phenotypic markers and marker intensities for each cell.
- *regionsprops* folder contains files for each patient, with X and Y coordinates and phenotypes for each cell.
- NIPMAP takes SampleIDs as a list to iterate through images - therefore the files needed to be re-labelled with integers.
- the script ***preprocessing_files.py*** will rename columns as in the NIPMAP tutorial (https://github.com/jhausserlab/NIPMAP/blob/main/README.md#inputs) and create a main csv where names of the files are coded into integers (SampleIDs), as well as create the cellData.csv with marker intensities, where cell labels and SampleIDs are cross referenced. <br>

For each SampleID: <br>
- Y_centroid -> *y*, y coordinate of the cell in the image <br>
- X_centroid -> *x*, x coordinate of the cell in the image <br> 
- Object -> *label*, integer cell label (matches with labels in cellData.csv) <br>
- Phenotype -> *cell_type*, string that describes the cell type <br>

File names: TS-373_IMC_patientID_disease_ROI.csv -> patient<SampleID>_cell_positions.csv.

For cellData.csv:
- Object -> *cellLabelInImage*, integer cell label (matches with labels in patient<SampleID>_cell_positions.csv) 
- *SampleID* column is added based on the patient ID
- remaining columns - markers with respective intensities 


### Running the analysis 

Prior to anything, the GitHub repository of [NIPMAP](https://github.com/jhausserlab/NIPMAP/tree/main) was cloned to the workspace. 

Running: 

```bash
bash running_nipmap.sh
```

will submit a batch job to the cluster and execute **running_nipmap.py**, where the main spatial analysis is happening. 

Make sure enough data is allocated before running the script. 

### Parameters to choose: 
In the running_nipmap.py main script, several parameters had to be changed to account for the bone dataset. 

Parameters to change: 
- **CELLTYPES**: list of cell phenotypes present across all images 
- **ImageIDs**: list of all SampleIDs, created in the pre-processing step 
- **NBNICHES**: number of niche, 7 in current script, can be adapted
- **RADIUS**: radii of sites (calculated based on optimal coverage & amount of PC), 13µm is sufficient to capture spatial niches
- **NSITES**: number of sites sampled from images, 100 are enough to capture the variance in a 1000x1000µm image.
- **ROOT_DATA_PATH**: directory where processed csv files with region properties are stored (patient<SampleID>_cell_positions.csv)

### Working with nipmap 
The python script will generate these outputs as json files:
pca_sites.json" # pca object on sites elements
AA_sites.json" # archetype Analysis object based on sites cell abundance
ca_sites.json" # cell abundance of randomly generated sites
cells_niches.json" # sites centered on cells and niches weights

The R script *py_wrapper_nipmap.r* 

- paramters to change: 

- mention the function script 

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
