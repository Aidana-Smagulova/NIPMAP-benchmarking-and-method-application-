# README

This report summarises the results of my internship with the Schapiro lab where I helped to benchmark spatial tumor microenvironment mapping tools: [NaroNet](https://www.sciencedirect.com/science/article/pii/S1361841522000366) and [NIPMAP](https://www.nature.com/articles/s41467-023-42878-z). My contribution included running these tools using the code provided on the GitHub repositories and checking, whether they deliver the same results as reported in the publications. Moreover, I aimed at customising the tools to apply them on the [myeloma dataset](https://github.com/SchapiroLabor/myeloma_standal) from Lukas Hatscher, to extend the metadata with the spatial tumor microenvironment information. 

NaroNet application was not successful, as the conda environment instructions provided in the code repository were not sufficient; and the combination of both TensorFlow and PyTorch libraries made the implementation of the tool tricky. It was decided to instead focus on NIPMAP as a more approacheable and promising tool. 

Current repository serves as a guide to reproduce NIPMAP application on the myeloma dataset.

This document's contents include:
- intstruction and scripts to run NIPMAP on the myeloma images 
- comparison and interpretation of the method output & conclusions

This README file bases on the original [NIPMAP GitHub repository](https://github.com/jhausserlab/NIPMAP/tree/main) and [publication](https://www.nature.com/articles/s41467-023-42878-z) by El Marrahi et al. 2023

## Index 


## Method description 
NIPMAP - NIche Phenotype MAPping (NIPMAP) analysis from spatial multiplex data.

NIPMAP is a statistical method that works to discover phenotypic niches in spatial data using principles from the community ecology approach, PCA and archetypes analysis. 

For this method, each cell type is treated as a species in an ecological niche. In the ecological ecosystems there are different species occupying different niches. Similarly, there are different cell types occupying different histological niches. Therefore, using the ideas from community ecology field we can consider histological niches as clusters of cells sites with similar histological profile. Within each niche cell types have certain density, which is how abundant this cell type is per surface area. The whole image is sampled randomly and uniformly in small circular cites with a radius that is big enough to capture more than one cell and so that less principal components might be needed to discover variance in cellular compositione. Per image, 100 sites are taken. This way, covered area represents 30% of the whole image, which allows efficient computations and accurate niche identification.  


Fig.1


Abundancies of each cell types are calculated inside of these sampling cites. Based on this, an abundance matrix of cell types in sampling sites is created. This abundance matrix is then used to perform PCA (according to authors, 3 principal components are sufficient to capture 82% of variance in cellular composition in images with sampling site radius of 13µm). After that, the PCA space is fitted onto a simplex figure using archetypal analysis (AA). An archetype in this case is an extreme niche case, when only one cell profile is abundant in a niche to 100% (say, only cancerous or only immune cells). Archetypes amount is determined manually - but it was shown that 4 niches are enough to capture over 80% of variance in cellular composition between sampling sites.

On the simplex, every point is a sampling site, which can be represented by a weighted average of tips of the simplex = archetypes. Thus, sites on the endpoints of the simplex, will have a 100% of one archetype weight and 0% of all others. These are solid niche representations in the tissue. A site that lays a little further away, will be calculated from all 4 archetypes, and highest weight can be decisive for niche allocation.

Niches might occur several times on a tissue slide, so colocalisation of niches will result in interfaces - cellular sites that consist of cell profiles from both neighbouring niches. Interface cell composition can be calculated as the weighted average of two niches it consists of. In the PCA, if a site lays directly in between two endpoints on the simplex, it will be considered an interface - a region between two niches. 

Therefore, NIPMAP is said to be able to discover histological niches in multiplex images. 
 
## Requirements 

Tested on: OS: Red Hat Enterprise Linux 8.10 (Ootpa)

Use the [minimal_env.yml](https://github.com/Aidana-Smagulova/NIPMAP-benchmarking-and-method-application-/blob/main/minimal_env.yml) file to recreate a working conda environment for NIPMAP. 

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

Further description of all packages in the environment can be found in the [requirements.txt](https://github.com/Aidana-Smagulova/NIPMAP-benchmarking-and-method-application-/blob/main/requirements.txt) file. 

* R libraries:
    ```
    - R 4.1.3
    - RStudio
    Open RStudio and run: 
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

NIPMAP works on cell coordinates and metadata, rather than images themselves. Therefore, a list of csv files that include X and Y coordinates for each cell to caclulate Gaussian Kernel distribution is needed. Moreover, cell phenotypes are used to derive the cell type abundancies in sampling sites and take advantage of the community ecology approach. Cellular markers and their intensities are further needed to map the phenotypic niches to cell markers. 

Data preprocessing steps for the myeloma dataset: 

Needed: two folders (intensities & regionsprops) with csv files for each patient, sample and ROI. names in format: TS-373_IMC78_UB_001.csv (TS-373_IMC_patientID_disease_ROI.csv)
- *intensities* folder contains files for each patient ID with phenotypic markers and marker intensities for each cell.
- *regionsprops* folder contains files for each patient, with X and Y coordinates and phenotypes for each cell.
- NIPMAP takes SampleIDs as a list of integers to iterate through files - therefore the original patient IDs need to be swapped with integers and then saved in a summary file to be able to map encoded file names with actual file names in the future.

The script [**preprocessing_files.py**](https://github.com/Aidana-Smagulova/NIPMAP-benchmarking-and-method-application-/blob/main/python_scripts/processing_files.py) will preprocess original data files stored in the *updated_csvs* folder.
- the script will process *regionprops* folder first. It will rename columns for each individual patient file as in the NIPMAP [tutorial](https://github.com/jhausserlab/NIPMAP/blob/main/README.md#inputs) and encode their names into NIPMAP readable format. For this, single files will be assigned a number from 1 to 154 (total number of patients in the dataset) and renamed. File renaming: TS-373_IMC_patientID_disease_ROI.csv -> patient<SampleID>_cell_positions.csv.
- The script will also create a csv file **patient_file_summary.csv** where original names of the files and NIPMAP readable names are stored for future referencing.
- Moreover, running the python script will create the **cellData.csv** file with marker intensities for each patient (from the *intensities* folder), where cell labels and SampleIDs are cross referenced with patient files (from the *regionprops* folder). It does so by concatenating the patient files with marker intensities together, matching the patient ID with the Sample ID from the **patient_file_summary.csv** and puts them in the SampleID column alongside with original cell IDs. <br>

Change of column names for each patient (sample) file: <br>
- Y_centroid -> *y*, y coordinate of the cell in the image <br>
- X_centroid -> *x*, x coordinate of the cell in the image <br> 
- Object -> *label*, integer cell label (matches with labels in cellData.csv) <br>
- Phenotype -> *cell_type*, string that describes the cell type <br>

Change of column names for **cellData.csv** file:
- Object -> *cellLabelInImage*, integer cell label (matches with labels in patient<SampleID>_cell_positions.csv) 
- *SampleID* column is added based on the patient ID
- remaining columns - markers with respective intensities

For the NIPMAP output post-processing R script **py_wrapper_nipmap.R** a **proteins_by_frame.csv** needs to be created. 
- it is a csv with Frame (index), Biomarker (list of all phenotypic marker) and Purpose (either "Element", "Background", "Functional" or "Lineage") columns
- will be needed for plotting in R
- save under *./phenotypes_niches/data/proteins_by_frame.csv*

### Running the analysis in Python 

Prior to anything, the GitHub repository of [NIPMAP](https://github.com/jhausserlab/NIPMAP/tree/main) was cloned to the workspace. 

```bash
git clone https://github.com/jhausserlab/NIPMAP.git
```

NIPMAP includes running the main.py script in terminal and then taking the json outputs to R for plotting and visualising. 

Executing main python script  **running_nipmap.py** will run calculate cell type abundancies in sampling sites, generate an abundance matrix, run PCA, archetypal analysis, and calculate niche weights for sampling sites. 

Running the running_nipmap.sh from inside of the cloned NIPMAP repository:

```bash
bash running_nipmap.sh
```

will submit a batch job to the cluster and execute **running_nipmap.py**, where the main spatial analysis is happening. 

*Make sure enough memory is allocated before running the script. In my case, 1TB was enough to run NIPMAP on 154 samples. Depending on the amount of niches and sample files, it can take up to one hour.*

```
#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --mem=1000gb
#SBATCH --export=NONE
```

 **running_nipmap.py** output:<br>
pca_sites.json – a PCA object on all sampling sites<br>
AA_sites.json – Archetype Analysis object based on sites cell abundance<br>
ca_sites.json – file with cell abundancies of randomly generated sampling sites<br>
cells_niches.json – file with sites centered on cells and niches weights<br>
params.json – file with parameters used to run the script (cell types, ImagesIDs, radius size, number of niches, number of sampling sites, etc)<br>
sites_cells_archs.csv – a big dataframe with niche weight, SampleIDs, and phenotypic markers for each cell ID

In comparison to the original **main_nipmap.py** script, several parameters had to be changed to account for the myeloma dataset. 

- **CELLTYPES**: list of cell phenotypes present across all images (same cell type names as in the processed regionprops csv files)
- **ImageIDs**: list of all SampleIDs, created in the pre-processing step (same SampleIDs as in the processed regionprops csv files)
- **NBNICHES**: number of niches, 7 in current script, can be adapted
- **RADIUS**: radii of sampling sites = 13µm (chosen to be the best for 3 principal components for the PCA)
- **NSITES**: number of sites sampled from images = 100 (allows to uniformly sample the image)
- **ROOT_DATA_PATH**: directory where processed csv files with region properties are stored (TMENS_analysis/data/cell_positions_data/...)

Both python script to run NIPMAP and preprocess data are available in current repository in the python_scripts folder. 

## Post-processing and plotting in R 

JSON files are then fed into the R script **py_wrapper_nipmap.r** (after installing necessary packages) 

The R script will: 
- read the json files
- plot the PCA space with the sampling sites fitted onto a simplex
- write a csv with cell types densities in each of the archetypes (nicheComposition.csv) and plot it into a barplot (barplotNiches.pdf)
- use cellData.csv to map niches & phenotypic markers in the cellPhenotypesNiches.csv
- plot heatmaps with niches & phenotypic markers correlation
- allocates niches to each sampling site and writes to sites_niches_all_patients.csv
- writes niches_ct_density for each of the marker density in all 7 niches
- writes csv Niches_whichcells.csv with certain cell types enriched in niches (basically niche composition)

Parameters to change: 
- data input (json files) & output paths (where all the resulting csv's should be written)

eg: jsonparams <- fromJSON(file="./lukas_data/updated_output/params.json")

- phenotypic markers (a list) for the mapping

line 124 of the **py_wrapper_nipmap.r** script: 

MARKERS <- c('CD38', 'Perilipin', 'Vimentin', 'B4GALT1', 'MPO',
             'CathepsinK', 'ATP5A', 'RUNX2', 'HIF1A', 'CD11b', 'CD45', 'CS', 'CD11c',
             'CD36', 'CD4', 'CD34', 'CD68', 'IL32', 'IDO', 'CD8', 'GranzymeK',
             'PKM2', 'IRF4', 'GLUT1', 'GranzymeB', 'Ki67', 'CollagenTypeI', 'CD3',
             'CPT1A', 'CD98', 'HLA-DR', 'ST6GAL1', 'CD138')
- niche names in the NichesNames variable

NichesNames <- c("a1"="niche1","a2"="niche2","a3"="niche3","a4"="niche4")

- use the custom proteins_by_frame.csv

Another R script **functions_phenotypes_tmens.r** contains all necessary functions to produce the csv's and plots in the py_wrapper_nipmap.r script. In this one, several functions had to be debugged and customised to fit the myeloma dataset. Therefore, both scripts are fixed and provided in the current repository with updated parameters. Right now scripts are written to run NIPMAP on 7 niches, but it can be fixed with parameters above. 

The functions script needs to be saved under *./phenotypes_niches/functions_phenotypes_tmens.r* path. 

The main plotting R script can be saved in the same directory as the **running_nipmap.py** script. 

Custom R script **creating_csv_for_plots.R** will: 
- use the same function as the **py_wrapper_nipmap.r** to read JSON files (primarily the cells_niches.json with sampling sites and their corresponding niche weights)
- open the json file and lay it out as a tibble 
- iterate through ./TMENS_analysis/data/cell_positions_data directory, where all patient csv files with cell coordinates, SampleIDs and phenotypes are
- create a combined dataframe with the data from all patients 
- match the csv files' SampleIDs with the SampleIDs from the json files 
- left-join the combined dataframe with patient data with the json file, based on cell_id and SampleID column
- as a result we get one big csv where each cell is a row and columns contain the information on cell coordinates, phenotype, SampleID it belongs to, plus niche and interface weights for all 7 niches and 21 interfaces

The resulting **cells_niches_coordinates_interfaces.csv** is used to create custom plots using the Jupyter Notebooks. 

### Results / Interpretation 

- **nipmap_result_evaluation.ipynb** plots niche cell-type correaltion & compares 7 NIPMAP niches to 7 k-means clusters


THE INTERFACE ISSUE 
ASSIGNING NICHES TO CELLS 


El Marrahi, A., Lipreri, F., Kang, Z. et al. NIPMAP: niche-phenotype mapping of multiplex histology data by community ecology. Nat Commun 14, 7182 (2023). https://doi.org/10.1038/s41467-023-42878-z
