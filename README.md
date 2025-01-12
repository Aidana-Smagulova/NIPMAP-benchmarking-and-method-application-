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

NIPMAP is a spatial mathod that works to discover niches (clusters of cells with similar profile) using principles from the community ecology approach, PCA and archetypes analysis. 
For this method, each cell type is treated like a species in an ecological niche. The whole image is sampled randomly and uniformly in small circular cites (radius can be decided manually), and abundancies of each cell types are calculated inside of the sampling cites. Based on this, an abundance matrix of cell types in sampling sites is created. 
This abundance matrix is then used to perform PCA (paper shows that 3 components are enough to capture the variance). After that, the PCA space is fitted onto a simplex figure using archetypal analysis (AA). An archetype in this case is a niche extreme case, when only one cell profile is abundant in a niche to 100% (say, only cancerous or only immune cells). Therefore, simplex tips are archetypes. If a sampling site falls on a samplex tip, the archetype weight will be 100% for this exact site. However, if a sampling site falls on an edge of the figure, it will not be allocated to a single archetype with 100% certainly but rather will represent a mixture of two near-lying archetype weights in percent. Sampling sites right in between of two archetypes (50% for each of them) are determined to be interfaces - regions that lay between any of the two niches. 


Within each niche, each cell type has a specific density, defined as the abundance of cells of that type per surface area of the niche
The niche recurs over the tissue section so that a limited number of niches is sufficient to capture the tissue’s cellular architecture, by piecing niches together. 
 The histological niches of the tissue are then revealed as clusters of sites with similar cellular composition
 Sites will form clusters of cellular composition if histological niches occupy distinct areas of the tissue with few interfaces, so that sites belong to only one niche (Fig. 1b). But, when niches colocalize and form larger interfaces, many sites lie at the interface of niches. Because the cellular composition of these sites is a mix (a weighted average) of the corresponding niches, no clear clusters can be distinguished by scattering the cellular composition of sites

 Sites are points positioned on axes that represent cellular composition. Because there are 17 cell types, 17 axes (dimensions) are needed in principle. This creates a representation challenge as human intuition is limited to 3 dimensions. However, two principles decrease the number of axes required to interpret tissue architecture. First, the abundance of certain cell types varies little across sites and thus contributes little to tissue architecture so that the corresponding axes can be neglected. Second, the abundance of cell types can correlate across sites—for example, because the cells cooperate in performing a tissue function—so that these cell types can be grouped into a single axis. The axes that optimally capture site cellular composition can be determined automatically by PCA, following the community ecology approach.

To interpret tissue architecture, it is important to set the radius of sampling sites to an appropriate size. Sites need to be large enough to capture local coordination in cellular composition and small enough to avoid blurring this coordination across different niches. To determine an appropriate radius, local coordination was quantified by the number of axes—principal components (PCs)—needed to capture spatial variation in cellular composition. When sampling sites are too small, they cover only one cell at a time: there is little covariance in the cellular composition of sites, and many axes (PCs) are thus needed to capture the cellular composition of sites. Increasing the radius of sites to include neighboring cells reveals covariance structures so that a smaller number of axes (PCs) is sufficient in capturing site cellular composition.


When sites have a 25μm radius, three PCs are enough to capture 82% of the variance in site cellular composition 

 A site radius of 25μm implies that cellular coordination emerges at a length scale of 2–4 cells. Increasing the site radius beyond 25 μm uncovered little novel covariance. We thus set the sampling site radius to 25 μm.

 Scattering sites on three PCs revealed no clear clusters (Fig. 2c). Instead, sites described a continuum with the shape of a 3D simplex: a pyramid with a triangular basis. This observation has significance for interpreting tissue architecture. Any point within a simplex can be described as a weighted average of the endpoints that define the simplex (Fig. 2d). Thus, observing that sites are constrained by the geometry of a 3D simplex implies that local cellular composition of the tissue is a mix (weighted average) of four histological niches, the endpoints of the 3D simplex (Fig. 2d). Sites close to endpoints represent cores within the niches, whereas sites halfway between endpoints localize at the interface between two niches (Fig. 2d). This interpretation generalizes the two-niches-and-interface interpretation of continua of site composition introduced in Fig. 1c to more than two niches.

 We note that niches were determined by collecting 100 sites per sample, so that the total area covered by sites represents 30% of the image area. Such a sampling intensity is sufficient to accurately identify niches while speeding up computations

 We find that 4 community ecology-based niches capture 82% of the inter-site variance in cellular composition while 4 clustering-based niches capture 58% of the variance

 
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

*Make sure enough data is allocated before running the script.*

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
