# README

Current repository serves as a guide to reproduce NIPMAP application on a bone disease dataset for benchmarking.

This document's contents include:
- intstruction and scripts to run NIPMAP on the bone cancer dataset (from Lukas Hatscher) 
- comparison and interpretation of the method output & conclusions

This README file bases on the original [NIPMAP GitHub repository](https://github.com/jhausserlab/NIPMAP/tree/main) and [publication](https://www.nature.com/articles/s41467-023-42878-z) by El Marrahi et al. 2023

## Index 


## Method description 
NIPMAP - NIche Phenotype MAPping (NIPMAP) analysis from spatial multiplex data.

Tissue microenvironments play a significant role in pathology. 

NIPMAP is a spatial mathod that works to discover niches using principles from the community ecology approach, PCA and archetypes analysis. 

For this method, each cell type is treated like a species in an ecological niche. Just like in ecological ecosystems there are different species occupying different niches, there are different cell types in histological niches. So, using the ideas from community ecology field, we can view histological niches as clusters of cells sites with similar histological profile. Within each niche, cell types have certain density, which is how abundant this cell type is per surface area. The whole image is sampled randomly and uniformly in small circular cites with a radius that is big enough to capture more than one cell and so that less principal components might be needed to discover variance in cellular compositione. Per image, 100 sites are taken. This way, covered area represents 30% of the whole image, which allows efficient computations and accurate niche identification.  

Abundancies of each cell types are calculated inside of the sampling cites. Based on this, an abundance matrix of cell types in sampling sites is created. This abundance matrix is then used to perform PCA (with image site radius of 13µm 82% of variance in cellular composition can be discovered using just 3 principal components). After that, the PCA space is fitted onto a simplex figure using archetypal analysis (AA). An archetype in this case is an extreme niche case, when only one cell profile is abundant in a niche to 100% (say, only cancerous or only immune cells). Archetypes amount is determined manually - but it was shown that 4 niches are enough to capture over 80% of variance in cellular composition between sampling sites.

On the simplex, every point is a sampling site, which can be represented by a weighted average of tips of the simplex = archetypes. Thus, sites on the endpoints of the simplex, will have a 100% of one archetype weight and 0% of all others. These are solid niche representations in the tissue. A site that lays a little further away, will be calculated from all 4 archetypes, and highest weight can be decisive for niche allocation.

Because niches might occur several times on a tissue slide, colocalisation of niches will result in interfaces - cellular sites with composition which is a weighter average of the neighboring niches. In the PCA, if a site lays directly in between two endpoints on the simplex, it will be considered an interface - a region between two niches. 

Therefore, NIPMAP is able to discover histological niches in multiplex images. 
 
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
- NIPMAP takes SampleIDs as a list of integers to iterate through files - therefore the original patient IDs need to be swapped with integers.

The script **preprocessing_files.py** will preprocess original data files stored in the *updated_csvs* folder.
- the script will process *regionprops* folder first. It will rename columns for each individual patient file as in the NIPMAP [tutorial](https://github.com/jhausserlab/NIPMAP/blob/main/README.md#inputs) and encode their names into NIPMAP readable format. For this, single files will be assigned a number from 1 to 154 (total number of patients in the dataset) and renamed. File renaming: TS-373_IMC_patientID_disease_ROI.csv -> patient<SampleID>_cell_positions.csv.
- The script will also create a csv file **patient_file_summary.csv** where original names of the files and NIPMAP readable names are stored for future referencing.
- Moreover, running the python script will create the **cellData.csv** file with marker intensities for each patient (from the *intensities* folder), where cell labels and SampleIDs are cross referenced with patient files (from the *regionprops* folder). It does so by concatenating the patient files with marker intensities together, matching the patient ID with the Sample ID from the **patient_file_summary.csv** and puts them in the SampleID column alongside with original cell IDs. <br>

Change of column names for each patient (sample) file: <br>
- Y_centroid -> *y*, y coordinate of the cell in the image <br>
- X_centroid -> *x*, x coordinate of the cell in the image <br> 
- Object -> *label*, integer cell label (matches with labels in cellData.csv) <br>
- Phenotype -> *cell_type*, string that describes the cell type <br>

Change of column names for cellData.csv file:
- Object -> *cellLabelInImage*, integer cell label (matches with labels in patient<SampleID>_cell_positions.csv) 
- *SampleID* column is added based on the patient ID
- remaining columns - markers with respective intensities

- mention that ./phenotypes_niches/data/proteins_by_frame.csv
- csv with Frame (index), Biomarker (list of all phenotypic marker) and Purpose (either "Element", "Background", "Functional" or "Lineage") columns
- needs to be created, will be needed for plotting in R

### Running the analysis 

Prior to anything, the GitHub repository of [NIPMAP](https://github.com/jhausserlab/NIPMAP/tree/main) was cloned to the workspace. 

```bash
git clone https://github.com/jhausserlab/NIPMAP.git
```

Running: 

```bash
bash running_nipmap.sh
```

will submit a batch job to the cluster and execute **running_nipmap.py**, where the main spatial analysis is happening. 

*Make sure enough memory is allocated before running the script.*

### Working with nipmap 

NIPMAP consists of running the main .py script and then taking the outputs to R for plotting and visualising. 

Executing main python script  **running_nipmap.py** will run calculate cell type abundancies in sampling sites, generate an abundance matrix, run PCA, archetypal analysis, and calculate niche weights for sampling sites. 

 **running_nipmap.py** output:<br>
pca_sites.json # pca object on sampling sites<br>
AA_sites.json # Archetype Analysis object based on sites cell abundance<br>
ca_sites.json # cell abundance of randomly generated sites<br>
cells_niches.json # sites centered on cells and niches weights<br>
params.json # file with parameters used to run the script (cell types, ImagesIDs, radius size, number of niches, number of sampling sites, etc)
sites_cells_archs.csv #niche weight, SampleIDs, and phenotypic markers for each cell ID

In comparison to the original **main_nipmap.py** script, several parameters had to be changed to account for the bone dataset. 

- **CELLTYPES**: list of cell phenotypes present across all images 
- **ImageIDs**: list of all SampleIDs, created in the pre-processing step 
- **NBNICHES**: number of niches, 7 in current script, can be adapted
- **RADIUS**: radii of sampling sites = 13µm
- **NSITES**: number of sites sampled from images = 100 
- **ROOT_DATA_PATH**: directory where processed csv files with region properties are stored (TMENS_analysis/data/cell_positions_data/...)

JSON files are then fed into the R script **py_wrapper_nipmap.r**

The R script will: 
- read json files
- plot the PCA space with the sampling sites fitted onto a simplex
- write a csv with cell types densities in each of the archetypes (nicheComposition.csv) and plot it into a barplot (barplotNiches.pdf)
- use cellData.csv to map niches & phenotypic markers in the cellPhenotypesNiches.csv
- plot heatmaps with niches & phenotypic markers correlation
- allocates niches to each sampling site and writes to sites_niches_all_patients.csv
- writes niches_ct_density for each of the marker density in all 7 niches
- writes csv Niches_whichcells.csv with certain cell types enriched in niches (basically niche composition)

Paramters to change: 
- data & output paths
- phenotypic markers (a list) for the mapping
- niche names in the NichesNames variable

functions_phenotypes_tmens.r contains all necessary functions to produce the csv's and plots in the py_wrapper_nipmap.r script. 

+ custom R script


i just took one json file to turn it into a giant csv and work further with it 



### results / interpretation / scatterplots 
COMPARING WITH K MEANS CLUSTERING 
THE INTERFACE ISSUE 
ASSIGNING NICHES TO CELLS 


El Marrahi, A., Lipreri, F., Kang, Z. et al. NIPMAP: niche-phenotype mapping of multiplex histology data by community ecology. Nat Commun 14, 7182 (2023). https://doi.org/10.1038/s41467-023-42878-z
