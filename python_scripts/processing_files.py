import pandas as pd
import os
import glob 

#initialise a list to store file names for mapping and a counter to iterate through the files 
patient_file_data = []
counter = 1 

#important! that files need to be ordered in the ascending manner from the lowest patient ID to the highest 

#specify output folders 
files_folder = './regionprops'
files_output_folder = './TMENS_analysis/data/cell_positions_data/'

os.makedirs(files_output_folder, exist_ok=True)

all_files = glob.glob(os.path.join(files_folder, "*.csv"))

#iterate through "regionprops" folder files to preprocess 
for filename in sorted(all_files):
    
    base_name = os.path.basename(filename)

    print(f"Processing file: {base_name}") #log 

    reg = pd.read_csv(filename)

    #renaming columns for compatibility with NIPMAP 
    reg.rename(columns={'Y_centroid': 'y', 'X_centroid': 'x', 'Object': 'label', 'Phenotype': 'cell_type', 'area': 'cell_size'}, inplace=True)

    #some columns are unnecessary, can be removed 
    reg = reg.drop(columns=['axis_major_length', 'axis_minor_length', 'eccentricity', 'distance_to_bone'])

    reg['SampleID'] = counter
    print(f"SampleID set to {counter} for {filename}")

    #columns are reordered
    new_column_order = ['x', 'y', 'label', 'cell_type', 'SampleID', 'cell_size']
    reg = reg[new_column_order]

    #counter is used to generate unique file names for each files 
    new_filename = f"patient{counter}_cell_positions.csv"

    #the OG file names and the generated ones are saved in a separate csv for further mapping 
    patient_file_data.append({
        "original_file_name": base_name,
        "patient_file_name": new_filename,
        "counter": counter
    })

    #writing new csv files in the output folder 
    reg.to_csv(os.path.join(files_output_folder, new_filename), index=False)

    #move on to the next file 
    counter += 1

print(f"Finished processing files")

#after iterating through all files and renaming them, save the overview in a separate csv 
summary_df = pd.DataFrame(patient_file_data)
summary_df.to_csv(os.path.join(files_output_folder, "patient_file_summary.csv"), index=False)

print("Patient file summary saved as patient_file_summary.csv")


#load the freshly written csv, as we're going to match the names with csv files from the "intensities" folder 
patient_file_data = pd.read_csv("TMENS_analysis/data/cell_positions_data/patient_file_summary.csv")

#initialise a list where files will be stored later - this is needed to concatenate all of the files into cellData.csv 
all_data = []

#specify data input & output directories 
intensities_directory = 'updated_csvs/intensities/' 

intensities_output_folder = "TMENS_analysis/data/cell_intensities_data/"
os.makedirs(intensities_output_folder, exist_ok=True)

#iterate through each file in the intensities directory
for filename in os.listdir(intensities_directory):
    if filename.endswith(".csv"):
        # Construct the full file path
        file_path = os.path.join(intensities_directory, filename)
        
        #check if the file exists also in the "regionprops" folder and match the OG name to then allocate the same encoded name to the same sample file in the intensity directory 
        match = patient_file_data[patient_file_data['original_file_name'] == filename]
        
        if not match.empty:
            sample_id = match.iloc[0]['counter']
            
            #read the intensity file 
            intensities = pd.read_csv(file_path)
            
            #add the matching SampleID column to the csv file 
            intensities['SampleID'] = sample_id
            
            #rename column in the file 
            intensities.rename(columns={'Object': 'cellLabelInImage'}, inplace=True)

            #reorder columns for better observability (plus NIPMAP is going to iterate through marker names later) 
            new_column_order = ['SampleID', 'cellLabelInImage', 'CD38', 'Perilipin', 'Vimentin', 'B4GALT1', 'MPO',
       'CathepsinK', 'ATP5A', 'RUNX2', 'HIF1A', 'CD11b', 'CD45', 'CS', 'CD11c',
       'CD36', 'CD4', 'CD34', 'CD68', 'IL32', 'IDO', 'CD8', 'GranzymeK',
       'PKM2', 'IRF4', 'GLUT1', 'GranzymeB', 'Ki67', 'CollagenTypeI', 'CD3',
       'CPT1A', 'CD98', 'HLA-DR', 'ST6GAL1', 'CD138']
        
            #use only available columns 
            intensities = intensities[[col for col in new_column_order if col in intensities.columns]]
            
            #add the file to a giant list of dataframes to merge further 
            all_data.append(intensities)
            
            print(f"Processed file: {filename} with SampleID: {sample_id}")
        else:
            print(f"No match found for file: {filename}")


#concatenate the dataframes with intensity data from the all_data list into one big dataframe 
merged_data = pd.concat(all_data, ignore_index=True)

#save the merged file as 'cellData.csv'
merged_data.to_csv(os.path.join(intensities_output_folder, 'cellData.csv'), index=False)
print(f"All CSV files have been merged and saved to {os.path.join(intensities_output_folder, 'cellData.csv')}")
