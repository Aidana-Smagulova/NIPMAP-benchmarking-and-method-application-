pkgs <- c("tidyverse","ggplot2","ade4","factoextra","plotly","igraph","reshape2","ggrepel","viridis","fdrtool","pheatmap","cluster","broom","pROC","ggpubr","devtools","ggridges")
install.packages(pkgs)
#install.packages("jsonlite")

gc()
rm(list=ls())
#.libPaths("/scratch/anissa.el/R_old/x86_64-redhat-linux-gnu-library/4.0")
.libPaths("/home/common/R")
library(rjson)
library(tidyverse)
library(purrr)

### SET WORKING DIRECTORY
dirName <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dirName)
source("./phenotypes_niches/functions_phenotypes_tmens.r")

jsonparams <- fromJSON(file="./7n_output/params.json")
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


file1 = "./7n_output/pca_sites.json" # pca object on sites elements
file2 = "./7n_output/AA_sites.json" # archetype Analysis object based on sites cell abundance
file3 = "./7n_output/ca_sites.json" # cell abundance of randomly generated sites
file4 = "./7n_output/cells_niches.json" # sites centered on cells and niches weights

#######---- Open .json files ----#######
json_data <- fromJSON(file=file1)
json_data2 <- fromJSON(file=file2)
json_data3 <- fromJSON(file=file3)
json_data4 <- fromJSON(file=file4)

##### LOAD OUTPUT OBJECTS
## Cell abundance in sites
sitesCellAb <- as_tibble(lapply(json_data3$cellAbSites,unlist))
write_csv(sitesCellAb,"./7n_output/sitesCA.csv")


niches <- paste0("a",as.vector(seq(1,NBNICHES,1)))
names(COLNICHES) <- niches
colNiches.hex <-unlist(lapply(COLNICHES, function(x){rgb(x[1],x[2],x[3],maxColorValue = 255)}))

## Archetypes coordinates in reduced PC space
Archs_3D <- do.call(cbind,lapply(json_data2$archs_coord,unlist))
## Projection of sites cell abundance in reduced PC space
pca3D <- matrix(unlist(json_data$PC_proj),nrow=15)[1:3,] 
#print(length(CELLTYPES))
plotly::plot_ly(x=pca3D[1,],
                y=pca3D[2,],
                z=pca3D[3,],
                type = "scatter3d", mode = "markers",
                marker = list(symbol = "triangle",size = 4),
                name="sites",
                mode = "text")%>%
  plotly::add_trace(x = Archs_3D[,1],
            y = Archs_3D[,2],
            z =Archs_3D[,3],
            type = "scatter3d",
            mode = "markers+text",
            text = niches,
            #textposition = c('top right','bottom right','top left','top right'),
            textposition = "top right",
            textfont = list(color = '#000000', size = 15),
            showlegend = TRUE,
            name = "niches",
            marker = list(color=~colNiches.hex,symbol = "star-diamond",size = 12),
            inherit = FALSE)%>%
  plotly::layout(scene = list(xaxis = list(title = "PC1"),
                       yaxis = list(title = "PC2"),
                       zaxis = list(title = "PC3")))


# arche_list <- c()
# for (i in seq(1,NBNICHES)) {
#   arche_list <- arche_list %>% append(paste("arch",i) %>% str_replace_all(" ",""))
# }
# 
# #######---- The pc composition of all niches ----#######
# 
# niche_pccompo <- read.csv(paste(RESULT_DATA_PATH,"niche_pcabundance.csv") %>% str_replace_all(" ","")) %>% 
#   dplyr::select(-c(X)) %>% 
#   cbind(arche_list = arche_list) %>% 
#   column_to_rownames(var = "arche_list")
# niche_alfas <- matrix(unlist(file2$alfas), ncol=length(arche_list)) %>% as.data.frame()
# colnames(niche_alfas) <- arche_list
# 
# plist <- Faceprojection(arche_list=, niche_pccompo, niche_alfas)
# plot <- do.call("grid.arrange", c(plist, ncol = 4))
# ggsave("./face_projection_niches.pdf", plot, height=10,width=12)

######--- NICHE IDENTIFICATION 

NichesCellProf <- do.call(cbind,lapply(json_data2$nichesCA,unlist))
rownames(NichesCellProf) <- CELLTYPES
colnames(NichesCellProf) <- niches
NichesCellProp <- NichesCellProf%>%t%>%as_tibble(rownames = NA)%>%
  rownames_to_column(var="archetype")%>%
  pivot_longer(cols=all_of(CELLTYPES),names_to="cell_type",values_to = "cell_density")
NichesCellProp[NichesCellProp<0] <-0
write_csv(NichesCellProp,"./7n_output/niche_cell_profile.csv")

barplot1 <- ggplot(data = NichesCellProp, aes(x = cell_type, y = cell_density,fill = archetype)) +
  geom_bar(stat = "identity",position = position_dodge(),width = 0.6) +
  scale_fill_manual(values = colNiches.hex)+
  theme(axis.text.x = element_text(angle = 90, vjust = .2)) + xlab ("") + ylab("cell density")
ggsave("./7n_output/barplot_niches.pdf",barplot1,height=3,width=4)
 

##########--------- NICHE-PHENOTYPE MAPPING --------##########
## Niches weights(proportions) of all cells from all images
niches<- paste0("a",as.vector(seq(1,NBNICHES,1)))
NINTERFACES <- 2
MARKERS <- c('CD38', 'Perilipin', 'Vimentin', 'B4GALT1', 'MPO',
             'CathepsinK', 'ATP5A', 'RUNX2', 'HIF1A', 'CD11b', 'CD45', 'CS', 'CD11c',
             'CD36', 'CD4', 'CD34', 'CD68', 'IL32', 'IDO', 'CD8', 'GranzymeK',
             'PKM2', 'IRF4', 'GLUT1', 'GranzymeB', 'Ki67', 'CollagenTypeI', 'CD3',
             'CPT1A', 'CD98', 'HLA-DR', 'ST6GAL1', 'CD138')

markers <- read.csv("./phenotypes_niches/data/proteins_by_frame.csv")%>%filter(Purpose=="functional")%>%pull(Biomarker)
cellsNiches <- as_tibble(lapply(json_data4$cells_niches,unlist))%>%
  rename_at(vars(matches("[0-9]")),~niches)%>%
  mutate(cell_id=as.numeric(cell_id))
#Make combinations of niches interfaces of order nIntf
getInterNiches <- function(nIntf,nbNiches){
  interfaces <- combn(paste0("a",as.vector(seq(1,nbNiches,1))),nIntf)
  coreIntf <- apply(interfaces,2,function(x) paste0(x,collapse=""))
  return(coreIntf)
}
coreIntf2 <- append(paste0("a",as.vector(seq(1,NBNICHES,1))),getInterNiches(NINTERFACES,NBNICHES))

# Get markers expression and niche weigths of cells
# df with following columns: cell_type, SampleID, cell_id, a1....an & interfaces, marker,value
cellsPhen.niches <- read.csv("./TMENS_analysis/data/cellData.csv",check.names=FALSE,header = TRUE, sep =',')%>%
  #dplyr::select(-c(cellSize,C,Na,Si,P,Ca,Fe,immuneCluster,Ta,Au))%>%
  #mutate(immuneGroup = recode(immuneGroup,`0`= 'None',`1`='Tregs', `2`='CD4-T',
                              #`3`='CD8-T', `4`='CD3-T', `5`='NK',
                             # `6`='B', `7`='Neutrophils', `8`='Macrophages', `9`='DC',
                             # `10`='DC / Mono', `11`='Mono / Neu', `12`='Other immune')) %>%
  #mutate(Group = recode(Group,`1`='Unidentified', `2`='Immune',
                      # `3`='Endothelial', `4`='Mesenchymal-like',
                       # `5` = 'Tumor',
                       # `6` = 'Keratin-positive tumor'))%>%
  #mutate(cell_type = ifelse(Group == 'Immune', cell_type<- immuneGroup,cell_type <- Group))%>%
  #dplyr::select(-c(tumorYN,tumorCluster,Group,immuneGroup))%>%filter(cell_type!="Unidentified")%>%
  #dplyr::rename(patient_id = SampleID)%>%
  dplyr::rename(cell_id = cellLabelInImage)%>%
  #mutate(H3K9ac_H3K27me3 = H3K9ac/H3K27me3)%>%
  left_join(cellsNiches%>%filter(cell_type!="Unidentified"),.,by=c("SampleID","cell_id"))%>%
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
  mutate(a6a7 = a6*a7)%>%
  pivot_longer(cols=all_of(MARKERS), 
               names_to="marker",values_to="value")


write_csv(cellsPhen.niches,"./cellPhenotypesNiches.csv")


source("./phenotypes_niches/functions_phenotypes_tmens.r")
CM <- correlation_niches_CM(markersCells.niches=cellsPhen.niches,Markers=markers,corrMeth="spearman",coreIntf2=any_of(coreIntf2),1/100,0.3,nbNiches=NBNICHES)

CM_with_rownames <- CM %>%
  rownames_to_column(var = "phenotypes/markers")

# Write to CSV
write.csv(CM_with_rownames, "./nichePhenMarker_correlation.csv", row.names = FALSE)

############# PLOT HEATMAPS AND TABLE OF SPATIAL CELL PHENOTYPES 
plot_heatmap_CT(CM.mat=CM,nichesIntf=any_of(coreIntf2),paste0(pathFigs,"/CMbyCells2.pdf"))
plot_heatmap_markers(CM.mat=CM,nichesIntf=any_of(coreIntf2),paste0(pathFigs,"/CMbyMarkers.pdf"))


##TABLE OF NICHE-ASSOCIATED CELL PHENOTYPES
#Archetype weights of randomly generated sites over TNBC patients images

#this is archetypes + weights on sites -> normalised to 1. 
archetypes_sites <- as.data.frame(do.call(cbind,lapply(json_data2$alfas,unlist)))

NichesNames <- c("a1"="niche1","a2"="niche2","a3"="niche3","a4"="niche4", "a5"="niche5", "a6"="niche6", "a7"="niche7")

#nice tibble where you can see that certain niches have more of a certain cell type
nichesCellAb <- NichesCellProf%>%t%>%as_tibble(rownames=NA)%>%rownames_to_column(var="archetype")%>%
  mutate(archetype=str_replace_all(archetype,NichesNames))%>%#str_replace_all(archetype,NichesNames))%>%
  column_to_rownames(var="archetype")%>%t()%>%
  as_tibble(rownames=NA)%>%rownames_to_column(var="cell_type")
# #Sort cell types per niche by decreasing cell abundance
#105 rows - 15 cell types x 7 niches -> so that we can see densities of each cell type in each niche 
nichesCA.sort <- nichesCellAb%>%
  pivot_longer(cols = as.vector(NichesNames),names_to="niche",values_to="cell_density")%>%
  group_by(niche)%>%arrange(desc(cell_density))

write_csv(as.data.frame(nichesCA.sort), "./niches_ct_density.csv")
colnames(archetypes_sites) <- niches

#sites with niche weights for all cell types in all patients (not pooled anymore)
archsSitesCellAb <- cbind(sitesCellAb,archetypes_sites)%>%
  pivot_longer(cols=all_of(CELLTYPES), names_to="cell_type", values_to="cell_density")

write_csv(as.data.frame(archsSitesCellAb), "./sites_niches_all_patients.csv")

### Get cell types enriched in each niche (take top 1% sites closest to niche)
archs.CT <- get_CT_enriched_all_archs(archsSitesCellAb,NichesNames)%>%#(archetypes_sites%>%t%>%as_tibble(rownames=NA),cellAbSites,thresh = 0.99)%>%
  group_by(niche)%>%
  #filter(!(cell_type %in% c("DC / Mono","CD3-T", "Mono / Neu", "Other immune","Unidentified")))%>%
  mutate(cell_type = paste(unique(cell_type),collapse="\n"))%>%
  distinct(niche,cell_type)
write_csv(archs.CT, "./niches_phenotypes_composition.csv")

source("./phenotypes_niches/functions_phenotypes_tmens.r")
## Get table of niches/interfaces-associated cell phenotypes
TableNichesPhenotypes(CM = CM,NichesCT = archs.CT,Niches.names = NichesNames,nichesCA.sorted = nichesCA.sort,pathFigs = pathFigs)
