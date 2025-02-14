# README

This report summarises the results of my internship with the Schapiro lab where I helped to benchmark spatial tumor microenvironment mapping tools: NaroNet by [Jiménez-Sánchez et al. 2022](https://www.sciencedirect.com/science/article/pii/S1361841522000366) and NIPMAP by [El Marrahi et al. 2023](https://www.nature.com/articles/s41467-023-42878-z). My contribution included running the tools and customising them to apply on the [myeloma dataset](https://github.com/SchapiroLabor/myeloma_standal) from Lukas Hatscher to possibly extend the metadata with the spatial tumour microenvironment information. 

NaroNet application was not successful, as the documentation provided in the code repository was not sufficient to run the tool and the combination of both TensorFlow and PyTorch libraries made the implementation of the tool tricky. It was decided to instead focus on NIPMAP as a more approachable and promising tool. 

The current repository serves as a guide to reproducing the NIPMAP application on the myeloma dataset.

This document's contents include:
- instruction and scripts to run NIPMAP on the myeloma images 
- interpretation of the method output & conclusions

This README file is based on the original [NIPMAP GitHub repository](https://github.com/jhausserlab/NIPMAP/tree/main) and [publication](https://www.nature.com/articles/s41467-023-42878-z) by El Marrahi et al. 2023


## Method description 
NIPMAP - NIche Phenotype MAPping (NIPMAP) analysis from spatial multiplex data.

NIPMAP is a statistical method for discovering phenotypic niches in spatial data using principles from community ecology approach, PCA, and archetypes analysis. 

For this method, each cell type is treated as a species in an ecological niche. In the ecological ecosystems, different species occupying different niches. Similarly, there are different cell types occupying different histological niches. Therefore, using the ideas from the community ecology field, we can consider histological niches as clusters of cell sites with similar histological profiles. Within each niche, certain cell types have a certain density, which is the abundance of a cell type per surface area. 

![this image shows a table with phenotypic niches and possible cell types they might include](https://github.com/Aidana-Smagulova/NIPMAP-benchmarking-and-method-application-/blob/main/figures/niche_example.png) 

*Fig.1 [El Marrahi et al. 2023](https://www.nature.com/articles/s41467-023-42878-z) Fig.1(d); Niches may include several cell types, but they to be similar in profile.* 
 
The whole image is sampled randomly and uniformly in small circular sites with a radius that is big enough to capture more than one cell and so that fewer principal components might be needed to discover variance in cellular composition. In the original publication, 100 sites were taken per image. This way, covered area represents 30% of a whole image, which allows efficient computations and accurate niche identification.  

![This image shows a phenotypic niche in grey laid over some cells in a spatial image](https://github.com/Aidana-Smagulova/NIPMAP-benchmarking-and-method-application-/blob/main/figures/niche_description.png)
*Fig.2 [El Marrahi et al. 2023](https://www.nature.com/articles/s41467-023-42878-z) Fig.1(a); Illustration of how niches consist of similar cell types that each have specific density. Density is cell type abundance per surface area unit.* 

Abundancies of each cell type are calculated inside of these sampling sites. Based on this, an abundance matrix of cell types in sampling sites is created. This abundance matrix is then used to perform PCA (according to authors, 3 principal components are sufficient to capture 82% of variance in cellular composition in images with a sampling site radius of 13µm). After that, the PCA space is fitted onto a simplex figure using archetypal analysis (AA). An archetype is an extreme niche case, when a niche consists of only one of the cell profiles (only cancerous or only immune cells). Archetypes number is determined manually - but it was shown that 4 niches are enough to capture over 80% of variance in cellular composition between sampling sites.

On the PCA space, every point is a sampling site. Each point on this simplex can be represented by a weighted average of tips of the simplex (archetypes).  Tips of the simplex are sampling sites that belong to only one niche to 100%, the weight of other niches in this site is 0. Thus, sites closer to one endpoint of the figure will have more of that archetype weight. However, a site that lays further away from the endpoints will be calculated from all 4 archetypes, resulting in the data-frame like format, where each cell has 4 weights for each niche. The highest niche weight for an individual cell can be decisive for niche allocation to this cell. 

![this image shows a PCA space, where each point is a sampling site and they are all situated on a simplex figure, where tips are archetypes](https://github.com/Aidana-Smagulova/NIPMAP-benchmarking-and-method-application-/blob/main/figures/pca.png)

*Fig.3 Representation of archetypal analysis result after fitting PCA space of sampling sites onto 7 archetypes. In the 3-dimensional space the PCA space takes the form of a simplex, and some archetypes overlap. PCA is generated from running NIPMAP on 154 samples of the myeloma dataset.*

Niches might occur several times on a tissue slide, so colocalisation of niches will result in interfaces - cellular sites that consist of cell profiles from both neighbouring niches. Interface cell composition can be calculated from the weighted average of two niches it consists of. In the PCA, if a site lays directly in between two endpoints on the simplex, it will be considered an interface - a region between two niches. 

![interface](https://github.com/Aidana-Smagulova/NIPMAP-benchmarking-and-method-application-/blob/main/figures/interfaces.png)
*Fig.4 [El Marrahi et al. 2023](https://www.nature.com/articles/s41467-023-42878-z) Fig.1 (d); This image represents how interfaces (regions between two niches) arise when niches colocalise and some cells fall of spaces between them, rather than falling onto one of the niches.* 

Therefore, NIPMAP is to be able to discover histological niches in multiplex images. 
 
## Requirements 

Tested on: OS: Red Hat Enterprise Linux 8.10 (Ootpa), Python ver. 3.7.13

Use the [minimal_env.yml](https://github.com/Aidana-Smagulova/NIPMAP-benchmarking-and-method-application-/blob/main/minimal_env.yml) file to recreate a working conda environment for NIPMAP. 

```bash
conda env create --file=minimal_env.yaml
```

Or run: 

```bash
conda create --name NIPMAP python=3.7.13
conda activate NIPMAP
pip install matplotlib scipy pandas==1.3.5 numpy>=1.17.3 scikit-learn seaborn qpsolvers==1.9.0
```

to recreate a working conda environment without a .yml file. 

Further description of all packages in the environment can be found in the [requirements.txt](https://github.com/Aidana-Smagulova/NIPMAP-benchmarking-and-method-application-/blob/main/requirements.txt) file, but it's not essential to create the NIPMAP conda environment. 

R: 
- R 4.1.3
- RStudio

Open RStudio and run
```R
pkgs <- c("tidyverse","ggplot2","ade4","factoextra","plotly","igraph","reshape2","ggrepel","viridis","fdrtool","pheatmap","cluster","broom","pROC","ggpubr","devtools","ggridges")

install.packages(pkgs)
```

## Method application

### Data: 

This is a description of how files should be ordered in the NIPMAP working directory. The TMENS_analysis folder should include the sample images csv files in the data/cell_positions_data directory; and a data/cellData.csv. Further descriptions below. 

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

### Data preprocessing steps for the myeloma dataset: 

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

After cloning the original repository, custom python and R scripts from this repository can be downloaded to the NIPMAP working directory, in order to run the tool.

NIPMAP includes running the main_nipmap.py script in terminal and then taking the .json outputs to R for plotting and visualising. Executing the main_nipmap.py script will calculate cell type abundancies in sampling sites, generate an abundance matrix, run PCA, archetypal analysis, and calculate niche weights for sampling sites. In case of the myeloma dataset, the main_nipmap.py is swapped with **running_nipmap.py** with adjusted parameters.  

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

In comparison to the original main_nipmap.py script, several parameters had to be changed to account for the myeloma dataset. 

- **CELLTYPES**: list of cell phenotypes present across all images (same cell type names as in the processed regionprops csv files)
- **ImageIDs**: list of all SampleIDs, created in the pre-processing step (same SampleIDs as in the processed regionprops csv files)
- **NBNICHES**: number of niches (7 in current script) can be adapted
- **RADIUS**: radii of sampling sites = 13µm (chosen to be the best for 3 principal components for the PCA)
- **NSITES**: number of sites sampled from images = 100q (allows to uniformly sample the image)
- **ROOT_DATA_PATH**: directory where processed csv files with region properties are stored (TMENS_analysis/data/cell_positions_data/...)

Both python script to run NIPMAP and preprocess data are available in current repository in the python_scripts folder. 

## Post-processing and plotting in R 

JSON files are then fed into the R script **py_wrapper_nipmap.r** (after installing necessary packages) 

The R script will: 
- read the json files
- plot the PCA space with the sampling sites fitted onto a simplex
- write nicheComposition.csv with cell types densities in each of the archetypes and plot it into the barplotNiches.pdf
- read cellData.csv to map niches & phenotypic markers in the cellPhenotypesNiches.csv
- plot heatmaps with niches & phenotypic markers correlation
- write sites_niches_all_patients.csv with niche weights for all cell types in all patients 
- write niches_ct_density.csv with densities of all cell types in all 7 niches (pooled for all patients) 
- write csv Niches_whichcells.csv with only highly enriched cell types in 7 niches (basically niche composition)

Parameters to change: 
- data input (json files) & output paths (where the resulting csv's should be written)
  ```R
  jsonparams <- fromJSON(file="./lukas_data/updated_output/params.json")
  ```

- phenotypic markers (a list) for the mapping
  line 124 of the **py_wrapper_nipmap.r** script:
  ```R
  MARKERS <- c('CD38', 'Perilipin', 'Vimentin', 'B4GALT1', 'MPO',
             'CathepsinK', 'ATP5A', 'RUNX2', 'HIF1A', 'CD11b', 'CD45', 'CS', 'CD11c',
             'CD36', 'CD4', 'CD34', 'CD68', 'IL32', 'IDO', 'CD8', 'GranzymeK',
             'PKM2', 'IRF4', 'GLUT1', 'GranzymeB', 'Ki67', 'CollagenTypeI', 'CD3',
             'CPT1A', 'CD98', 'HLA-DR', 'ST6GAL1', 'CD138')
  ```

- niche names in the NichesNames variable
  ```R
  NichesNames <- c("a1"="niche1","a2"="niche2","a3"="niche3","a4"="niche4", "a5"="niche5", "a6"="niche6", "a7"="niche7")
  ```

- use the custom **proteins_by_frame.csv**

Another R script **functions_phenotypes_tmens.r** contains all necessary functions to produce the csv's and plots in the **py_wrapper_nipmap.r** script. 

In **functions_phenotypes_tmens.r**, several functions had to be debugged and customised to fit the myeloma dataset. Therefore, both scripts are fixed and provided in the current repository with updated parameters. Right now scripts are written to run NIPMAP on 7 niches, but it can be fixed with parameters above. 

The functions script needs to be saved under *./phenotypes_niches/functions_phenotypes_tmens.r* path. 

The main plotting R script can be saved in the same directory as the **running_nipmap.py** script. 

Custom R script **creating_csv_for_plots.R** will: 
- use the same function as the **py_wrapper_nipmap.r** to read json files (primarily the cells_niches.json with sampling sites and their corresponding niche weights)
- open the json file and lay it out as a tibble 
- iterate through *./TMENS_analysis/data/cell_positions_data* directory, where all patient csv files with cell coordinates, SampleIDs and phenotypes are
- create a combined dataframe with the data from all patients 
- match the csv files' SampleIDs with the SampleIDs from the json files 
- left-join the combined dataframe with patient data with the json file, based on cell_id and SampleID column
- as a result we get one big csv where each cell is a row and columns contain the information on cell_ids, cell coordinates, phenotype, SampleID it belongs to, plus niche and interface weights for all 7 niches and 21 interfaces

The resulting **cells_niches_coordinates_interfaces.csv** is used to create custom plots using the Jupyter Notebooks. 

Further analysis is performed with this format of csv: 

![csv table with cells and niches](https://github.com/Aidana-Smagulova/NIPMAP-benchmarking-and-method-application-/blob/main/figures/cell_niche_csv.png)

*Fig.5 Main NIPMAP output after post-processing in R: csv table written from the cells_niches.json file, where each column is an entry relevant to an individual cell, and rows contain information such as 7 niche weights, cell ID, Sample ID and cell type for each niche.* 

### Results / Interpretation 

**nipmap_result_evaluation.ipynb** plots niche cell-type correaltion & compares 7 NIPMAP niches to 7 k-means clusters

Part of post-processing of NIPMAP results was unpacking the json files and presenting the output in a csv file that summarises the main points of NIPMAP. Mainly, individual niche weights for each cells. In the result, we aimed to have a big summary of cells and niches they most likely belong to, as well as cell coordinates, SampleIDs and cell types. The csv output from **cells_niches_coordinates_interfaces.csv** stores all 7 niches for each of the cells and 21 interfaces. The interface weights are calculated from individual niches - e.g. the interface weight of "a1a7" interface is caclualted by multiplication of niche weights "a1" and "a7". 

![csv table with cells and niches these cells most likely belong to](https://github.com/Aidana-Smagulova/NIPMAP-benchmarking-and-method-application-/blob/main/figures/cell_niche_interface_csv.png)

*Fig.6 CSV table where each cell is allocated to a niche AND an interface it most likely belongs to. Cells are allocated to the niches/interfaces with highest weight using a rowmax.*

The final result investigation also includes looking into composition of each niche - which cell types are allocated to which niches. For this, the frequencies of each cell type can be calculated from the **cells_niches_coordinates_interfaces.csv** and plotted in a heatmap manner, where warmer and brighter color represents higher cell type density in the given niche. 

![7 NIPMAP niches vs cell types heatmap](https://github.com/Aidana-Smagulova/NIPMAP-benchmarking-and-method-application-/blob/main/figures/niche_heatmap.png) 

*Fig.7 Heatmap of 7 NIPMAP niches composition, after running it on 154 samples of the myeloma dataset, where each row is a cell type and the abundance in a niche (column).* 

From here it can be seen that NIPMAP find niches that are enriched in different cell types. NIPMAP found a tumor niche a1 with plasma cells and several immune niches including a Neutrophil, Macrophage, MPO+, and T-cell niches. Moreover, NIPMAP niches even include several cell type that are similar: for example, niche 6 consist of T-cells, which shows that NIPMAP is capable of co-assotiating certain niches with very similar cell profiles. In the myelome dataset, some cell types are more abundant than other, e.g. plasma cells, neutrophils, or the unknown phenotype, and their density is reflected on the heatmap in comparison to less abundant cell types. 

After running the basic analysis and looking into niche composition, **scatter_plotting.ipynb** can be used to plots cells colored by niche & interfaces using **cells_niches_coordinates_interfaces.csv**. 

### Discussion

NIPMAP was employed to find 7 niches in the 154 images of the myeloma dataset and produce output in form that we can further work with. 

To evaluate the niche accuracy, NIPMAP output was compared to 7 k-means clusters from scimap that were ran on the same dataset. 

The correlation heatmap (Fig.8 and Frig.9) suggests that certain niches and k-means cluster might overlap in their cell type compositions. Which helps to see that certain NIPMAP niches might correlate with k-means clusters.

![7 k-means clusters from scimap vs cell types heatmap](https://github.com/Aidana-Smagulova/NIPMAP-benchmarking-and-method-application-/blob/main/figures/kmean_cluster_heatmap.png) 

*Fig.8 Heatmap of 7 k-means cluster composition, after running it on 154 samples of the myeloma dataset, where each row is a cell type and the abundance in a niche (column).* 

![correlation heatmap of 7 NIPMAP niches vs 7 k-means clusters from scimap](https://github.com/Aidana-Smagulova/NIPMAP-benchmarking-and-method-application-/blob/main/figures/heatmap_correlation.png)

*Fig.9 Heatmap of 7 k-means clusters from scimap, ran on 154 samples of the myeloma dataset. Each row is a cell type and each entry is its abundance in a certain k-mean cluster (column).* 

After examining niche composition we wanted to look into how these niches co-localise spatially and whether we can find correlation between niches and their interfaces. The main goal was to show that interfaces do in fact emerge in areas between niches. However, out data had the problem of allocating one cell to either a niche or an interface. Considering the biology, one cell cannot belong to both at the same time. It would make sense to produce an output where each cell could be allocated to either a niche or an interface, but here we faced a definition issue: how do we decide, which cells can belong to interfaces and which to just niches? 

The simplest solution was to find a cut-off value for interfaces - a certain threshold of an interface weight after which it can be said that a cell belongs to an interval with high accuracy. Thus, the higher the interface value, the more likely a cell is to be allocated to this interval. Several cut-off values were suggested. As we supposed, the cell should be allocated to an interface, if two of the seven niches have 50% weights. Therefore, if the calculated the interface weight for this cell as the product niche weights, a maximum interface weight value should be 0.25, as interfaces calculated from other niches would have lower wights. Several cut-off points were thus considered 0.18, 0.2, 0.22, and 0.24. The problem was not solved this way, as the plotting showed, with high cut-off values there were very little cells left in some images that could catogorise as interfaces. What is more important, even cells with very highest interface weights were still not colocalising with niche borders and sometimes emerged in a middle of one niche, where there was no visible niche/niche interaction. 

Then it was decided to take 25% of the rows with highest values from the sorted interface weight distribution and plot in a scatterplot. This looked well in some samples but not the others, there were still unexpected interface cells emerging in the middle of niches. Therefore, there is currently no universal solution to the interface issue that would work for all samples equally well. 

Since interface weights are calculated from niche weights, NIPMAP's accuracy in niche discovery played a major in interface determination.  NIPMAP was tested on 154 images, 100 sampling spots taken from each of them with a radius of 13µm. Considering that images are 1000µx1000µ in size, the sampling intensity of each image is approximately 5%, while in the paper by El. Marrahi et al. 2023 the minimal sampling intensity NIPMAP was tested on is 30%. The radius size is appropriate to run a PCA with 3 principal components, so one possibility would be to increase the sampling site number from one image, thus increasing sampling intensity, and, maybe, discovering more reliable niches and more accurate niche weights for each cell. If the cells are allocated to niches more reliably, maybe there could be an interface cut-off value determined that could work for at least most of the samples.  

One of the possible explanations is that NIPMAP algorithm could not work on the myeloma dataset due to factors unique to these images, for example tissue heterogeneity, colocalisation of many different cell types in same spots, unphenotyped cells like bones, etc. It is possible. It is worth discussing whether a smaller benchmark with less samples and very simple tissue architechtures could work better in order to estimate the accuracy of the tool. 

In conclusion, NIPMAP is a running tool that combines statistics and community ecology approach to discover tumor microenvironments. Initially, it was easy to implement and required little parameter change in order to optimise it for the myeloma dataset.

The question of how to quantify the accuracy of niche/interface allocation other than taking the rowmax remains. It is of crucial importance to solve the interface issue in order to draw any conclusions about the tool and its accuracy. Further analysis is also needed to validate and verify NIPMAP results, e.g. through cohort comparison between myeloma disease and normal samples and niches that are being discovered there. Plus, there are some things yet to try, for example the [power analysis](https://github.com/jhausserlab/NIPMAP/tree/main/power_analysis) script in the NIPMAP GitHub repository to detect rare niches, run sampling intensity analysis for the myeloma dataset, or run NIPMAP with different parameters, such as different amount of principal componentes, site sampling, site radii and/or niche number. 

[El Marrahi, A., Lipreri, F., Kang, Z. et al. NIPMAP: niche-phenotype mapping of multiplex histology data by community ecology. Nat Commun 14, 7182 (2023). https://doi.org/10.1038/s41467-023-42878-z](https://www.nature.com/articles/s41467-023-42878-z)
