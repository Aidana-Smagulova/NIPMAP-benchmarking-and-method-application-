import os
import sys
myDir = os.getcwd()
#print(sys.path)
module_path = myDir + "/TMENS_analysis/"
if module_path not in sys.path:
    sys.path.append(module_path)
#print(sys.path)
import pandas as pd
import numpy as np
from src.utils.archetypes import ArchetypalAnalysis
from src.CellAbundance import CellAbundance, join_abundance_matrices, generate_abundance_matrix
from src.utils.equations import compute_cells_niches_weights,get_niches_cell_abund
from src.utils.visualization import plot_cells_positions,radius_pc_all_variance,plot_cells_positions,plot_cells_TI_border
from sklearn.decomposition import PCA
import json
import shutil
from shutil import make_archive



#### ARGUMENTS
# CELLTYPES = list(sys.argv[1].split(",")) #['CD8-T', 'Other immune', 'DC / Mono', 'CD3-T', 'B', 'NK', 'Keratin-positive tumor', 'Tumor','CD4-T', 'Mesenchymal-like', 'Macrophages', 'Endothelial', 'Tregs', 'Unidentified', 'DC', 'Mono / Neu','Neutrophils']
# ImageIDs = [int(i) for i in sys.argv[2].split(",")] #[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37,38, 39, 40, 41]
# NSITES = int(sys.argv[3]) #100
# RADIUS = int(sys.argv[4]) #25
# NBNICHES = int(sys.argv[5]) #4
# METHOD = int(sys.argv[6]) #"gaussian"
# XSIZE= int(sys.argv[7])
# YSIZE = int(sys.argv[8])
# ROOT_DATA_PATH = sys.argv[9]  #"./TMENS_analysis/data/cell_positions_data"
# ROOT_OUTPUT_PATH = sys.argv[10] #"./TMENS_analysis/output"
# COLARCHS = np.array([[255, 0, 223],[255,0,0],[70,203,236],[0,0,0]]) #TODO  sys.argv[9]
CELLTYPES = ['Unknown', 'Osteocyte', 'Plasma Cells/MM cells', 'Neutrophils',
 'Macrophages/Monocytes', 'MPO+', 'CD8+Tcells',
 'activated Macrophages/Monocytes', 'CD68+', 'HSCs', 'CD4+Tcells',
 'Osteoblasts', 'Dendritic Cells', 'Adipocytes', 'Endothelial cells']
ImageIDs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154]
NSITES = 100
RADIUS = 10
NBNICHES = 7
METHOD ="gaussian"
XSIZE= 1000
YSIZE = 1000
ROOT_DATA_PATH = "./TMENS_analysis/data/cell_positions_data"
ROOT_OUTPUT_PATH = "./7n_output"

COLARCHS = np.array([[255, 0, 223],[255,0,0],[70,203,236],[0,0,0], [0, 255, 0], [0, 0, 255], [255, 255, 0], [255, 128, 0]])

## create figs directory for niches
#path_figs = "./7n_output/figs"
#path_toFigs = os.path.join(myDir,path_figs)


if __name__ == "__main__":
  
  #####----- RADIUS ANALYSIS -----####
  rad = np.linspace(np.log(5), np.log(190), num=10)# get interval of radiuses
  radius = np.rint(np.power(np.e, rad).astype(float))
  METHOD = 'gaussian'
  
  expl_var_ratio_gauss = {}
  # Get explained variace of each PC on sites cell abundance interating on radius size of sites
  for r in radius:
      print(r)
      gaussian_count_list = generate_abundance_matrix(CELLTYPES, ImageIDs, NSITES, r, center_sites_cells=False,method=METHOD,image_x_size=XSIZE,image_y_size=YSIZE, snr=3, root=ROOT_DATA_PATH)
      sites, patient_ids,s_ids, _ = join_abundance_matrices(gaussian_count_list)
      pca = PCA()
      pc = pca.fit_transform(sites)
      expl_var_ratio_gauss[r] = np.cumsum(pca.explained_variance_ratio_)
  radius_pc_all_variance(expl_var_ratio_gauss,radius_lim=13,nPC_lim=3,cells_number=len(CELLTYPES)+1,save_fig=True, path_fig="./plot_rad_var_gauss.svg")
  RADIUS = 13
  #int(input('Enter the size (in micrometers) of the radius of the sites (should be int):\n'))
  #print(type(RADIUS))

  #####----- GENERATE SITES AND COMPUTE CELL ABUNDANCE ----#####
  CellAb_list = generate_abundance_matrix(CELLTYPES, ImageIDs, NSITES,RADIUS,method=METHOD,random_seed=1022,snr=3,center_sites_cells=False,root=ROOT_DATA_PATH)
  sites, patients_ids,sites_ids, _ = join_abundance_matrices(CellAb_list)
  CellAb_df = pd.DataFrame()
  print("Generating sites with cell abundance...")
  for ca in CellAb_list:
      abundance_df = pd.DataFrame(ca.abundance_matrix,columns = CELLTYPES)
      abundance_df['site_id'] = np.arange(len(abundance_df))
      abundance_df['patient_id'] = ca.patient_id
      CellAb_df = CellAb_df.append(abundance_df)
  CellAb_df = CellAb_df.reset_index()
  #####----- PCA ON SITES ABUNDANCE ----#####
  print("Dimension reduction of sites cell abundances...")
  pca_obj = PCA()
  pc_proj = pca_obj.fit_transform(sites)
  #####----- ARCHETYPE ANALYSIS ----#####
  print("Finding niches...")
  AA = ArchetypalAnalysis(n_archetypes = NBNICHES,
                     tolerance = 0.001,
                     max_iter = 200,
                     random_state = 0,
                     C = 0.0001,
                     initialize = 'random',
                     redundancy_try = 30)
  AA.fit_transform(pc_proj[:,:NBNICHES-1])
  print(str(NBNICHES)+" niches found!")
  # 
  niches_cell_profile = get_niches_cell_abund(sites,pca_obj,AA,NBNICHES-1)
  #pd.DataFrame(niches_cell_profile)
  
  #if os.path.isdir(path_toFigs)==False:
  #  os.mkdir(path_toFigs)
  #ImageId = [18]
  # ##########----- SAVE PLOTS IN ZIP DIRECTORY ----##########
  
  print("Segmenting images into niches...")
  for i in ImageIDs:
    GRANULARITY = 5
    cell_data = pd.read_csv(ROOT_DATA_PATH+"/patient{}_cell_positions.csv".format(i))
    fig= plot_cells_TI_border(cell_data, CELLTYPES, patientID=ImageIDs,
                              data_niches_path ='../../7n_output/archetypes_sites_centered_all_cells_gaussian.csv',
                              intf=["arch2","arch3"],normalize=True,quant=0,path_fig="../../7n_output/figs/cells_patient"+str(patientID)+"_InflCancer_bordersN.svg")
  
  # shutil.make_archive("/figs_niches","zip", path_toFigs)
  
  #####----- GENERATE SITES CENTERED ON CELLS AND THEIR NICHE WEIGHTS ----#####
  print("Computing cells' niche weights, the operation might take some time...")
  CellAbCC_list = generate_abundance_matrix(CELLTYPES, ImageIDs, NSITES,RADIUS,method=METHOD, snr=3,center_sites_cells=True, border=False,root=ROOT_DATA_PATH)
  sitesCC, patients_ids2,sites_ids2, _ = join_abundance_matrices(CellAbCC_list)
  CellAbCC_df = pd.DataFrame()
  for ca in CellAbCC_list:
    df_ca = ca.get_site_cell_id_df()
    df_ca['patient_id'] = int(ca.patient_id)
    CellAbCC_df = CellAbCC_df.append(df_ca)
  CellAbCC_df = CellAbCC_df.reset_index(drop = True)

  NichesProf = get_niches_cell_abund(sitesCellAb=sites,pcaSites=pca_obj,ArchObj=AA,nComp=NBNICHES-1)
  sites_alfa = compute_cells_niches_weights(niches=NichesProf,cellsSites=sitesCC,nbNiches=NBNICHES)
  sites_archs = pd.DataFrame(sites_alfa)
  sites_archs['SampleID'] = patients_ids2
  sites_archs["cell_id"] = sites_ids2[:,0]
  sites_archs["cell_type"] = sites_ids2[:,1]
  sites_archs["TOT_cell_dens"]= sitesCC.sum(axis=1)
  sites_archs.to_csv("./7n_output/sites_cells_archs.csv",index=False)

  #TODO add plot_cells_TI_border for all images and only for tumor-immune borders (+ one parameter)
  #TODO directory with figures of segmented images and zip it  

  ##########----- SAVE OUTPUTS IN CSV ----##########
  ## SAVE PCA AND Archetype Analysis OBJECTS
  print("Saving outputs in json files...")
  dict_pca = {"PC_proj":pc_proj.tolist(),"components":pca_obj.components_.tolist(),"expl_variance":pca_obj.explained_variance_.tolist(),
    "expl_var_ratio":pca_obj.explained_variance_ratio_.tolist(),"mean":pca_obj.mean_.tolist()}

  dict_AA = {"archs_coord": AA.archetypes.tolist(), "alfas": AA.alfa.tolist(),"nichesCA":niches_cell_profile.tolist()}

  dict_caSites = {"cellAbSites": CellAb_df.to_dict()}
  dict_caSitesCC = {"cells_niches": sites_archs.to_dict(),"cellAb_sitesCC": CellAbCC_df.to_dict()}
  dict_params = {"cellTypes":CELLTYPES,"ImageID":ImageIDs,"nbsites":[NSITES],"radiusSize":[RADIUS],"nbniches":[NBNICHES],"countMeth":[METHOD],"xsize":[XSIZE],"ysize":[YSIZE],"rootDataPath":[ROOT_DATA_PATH],"rootOutPath":[ROOT_OUTPUT_PATH],"colNiches":COLARCHS.tolist()}

  # Serializing json objects
  PCA_json = json.dumps(dict_pca, indent=4)
  AA_json = json.dumps(dict_AA, indent=4)
  caSites_json = json.dumps(dict_caSites,indent=4)
  cellsNiches_json = json.dumps(dict_caSitesCC,indent=4)
  params_json = json.dumps(dict_params,indent=4)

  # Writing to .json files
  with open("./7n_output/pca_sites.json", "w") as outfile:
    outfile.write(PCA_json)
  
  with open("./7n_output/AA_sites.json","w") as outfile2:
    outfile2.write(AA_json)
  
  with open("./7n_output/ca_sites.json","w") as outfile3:
    outfile3.write(caSites_json)

  with open("./7n_output/cells_niches.json","w") as outfile4:
    outfile4.write(cellsNiches_json)
    
  with open("./7n_output/params.json","w") as outfile5:
    outfile5.write(params_json)

