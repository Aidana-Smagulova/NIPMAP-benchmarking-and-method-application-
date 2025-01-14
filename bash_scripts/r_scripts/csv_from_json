pkgs <- c("tidyverse","ggplot2","ade4","factoextra","plotly","igraph","reshape2","ggrepel","viridis","fdrtool","pheatmap","cluster","broom","pROC","ggpubr","devtools","ggridges")
install.packages(pkgs)

.libPaths("/home/common/R")
library(rjson)
library(tidyverse)
library(purrr)

### SET WORKING DIRECTORY
dirName <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dirName)
source("./phenotypes_niches/functions_phenotypes_tmens.r")


#after transferring json files that the nipmap python script results in to the local directory, we can process them in R 
#here we read the json files 
jsonparams <- fromJSON(file="./lukas_data/7n_output/params.json")
CELLTYPES <-jsonparams$cellTypes
ImageIDs <- jsonparams$ImageID
NSITES <- jsonparams$nbsites
RADIUS <- jsonparams$radiusSize
NBNICHES <- jsonparams$nbniches
METHOD <-jsonparams$countMeth
W <-jsonparams$xsize
H <-jsonparams$ysize
ROOT_DATA_PATH <- jsonparams$rootDataPath
ROOT_OUTPUT_PATH <-jsonparams$rootOutPath
COLNICHES <- jsonparams$colNiches
pathFigs <- jsonparams$pathFigs


file1 = "./lukas_data/7n_output/pca_sites.json" # pca object on sites elements
file2 = "./lukas_data/7n_output/AA_sites.json" # archetype Analysis object based on sites cell abundance
file3 = "./lukas_data/7n_output/ca_sites.json" # cell abundance of randomly generated sites
file4 = "./lukas_data/7n_output/cells_niches.json" # sites centered on cells and niches weights

#######---- Open .json files ----#######
json_data <- fromJSON(file=file1)
json_data2 <- fromJSON(file=file2)
json_data3 <- fromJSON(file=file3)
json_data4 <- fromJSON(file=file4)


#define niches 
niches <- paste0("a",as.vector(seq(1,NBNICHES,1)))
names(COLNICHES) <- niches
colNiches.hex <-unlist(lapply(COLNICHES, function(x){rgb(x[1],x[2],x[3],maxColorValue = 255)}))

#create the CSV with 7 niches and 21 interface weights for each cell. 
plotting_niches <- as_tibble(lapply(json_data4$cells_niches,unlist))%>%
  rename_at(vars(matches("[0-9]")),~niches)%>%
  mutate(cell_id=as.numeric(cell_id))%>%
  filter(!(is.na(a1)| is.na(a2)| is.na(a3) | is.na(a4)))%>%
  mutate(a1a2 = a1*a2)%>%
  mutate(a1a3 = a1*a3)%>%
  mutate(a1a4 = a1*a4)%>%
  mutate(a1a5 = a1*a5)%>%
  mutate(a1a6 = a1*a6)%>%
  mutate(a1a7 = a1*a7)%>%
  mutate(a2a3 = a2*a3)%>%
  mutate(a2a4 = a2*a4)%>%
  mutate(a2a5 = a2*a5)%>%
  mutate(a2a6 = a2*a6)%>%
  mutate(a2a7 = a2*a7)%>%
  mutate(a3a4 = a3*a4)%>%
  mutate(a3a5 = a3*a5)%>%
  mutate(a3a6 = a3*a6)%>%
  mutate(a3a7 = a3*a7)%>%
  mutate(a4a5 = a4*a5)%>%
  mutate(a4a6 = a4*a6)%>%
  mutate(a4a7 = a4*a7)%>%
  mutate(a5a6 = a5*a6)%>%
  mutate(a5a7 = a5*a7)%>%
  mutate(a6a7 = a6*a7)



#define the directory containing the cell position data files
#this one contains processed CSV's with cell positions data in the format that nipmap understands (these stem from the updated anndata object)
directory <- "./lukas_data/full_dataset/cell_positions_data"

# iterate over each file in the directory
files <- list.files(directory, full.names = TRUE)

#combine all patient CSV files into one big file 
combined_data <- data.frame()


for (file in files) {  
  data <- read_csv(file, show_col_types = FALSE)
  
  #label column is renamed to cell_id
  data <- data %>% 
    rename(cell_id = label)
  
  #checking if sample id exists 
  if (!"SampleID" %in% colnames(data)) {
    stop(paste("SampleID column is missing in file:", file))
  }
  
  data <- data %>% 
    select(cell_id, SampleID, x, y)
  
  combined_data <- bind_rows(combined_data, data)
}

#here we left join the 4 columns celected from the cell position data files with the file containing niches & interfaces based on cell id and sample id 
#this way, each cell now has data about niche & interfaces weights as well as x and y coordinates (needed for plotting)
plotting_niches <- plotting_niches %>%
  left_join(combined_data, by = c("cell_id", "SampleID"))

write.csv(plotting_niches, "./lukas_data/7n_output/cells_niches_coordinates_interfaces_r10.csv")
