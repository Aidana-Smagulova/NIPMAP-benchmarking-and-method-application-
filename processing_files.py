import pandas as pd
import os
import glob 

# List to store patient file data
patient_file_data = []
counter = 1  # Initialize counter for unique IDs across all directories


# Output folders
files_folder = 'updated_csvs/regionprops'
files_output_folder = 'TMENS_analysis/data/cell_positions_data/'

os.makedirs(files_output_folder, exist_ok=True)

all_files = glob.glob(os.path.join(files_folder, "*.csv"))

# Process each file in the directory
for filename in sorted(all_files):
    
    base_name = os.path.basename(filename)

    print(f"Processing file: {base_name}")

        # Read the CSV file
    reg = pd.read_csv(filename)

        # Rename columns
    reg.rename(columns={'Y_centroid': 'y', 'X_centroid': 'x', 'Object': 'label', 'Phenotype': 'cell_type', 'area': 'cell_size'}, inplace=True)

        # Drop unnecessary columns
    reg = reg.drop(columns=['axis_major_length', 'axis_minor_length', 'eccentricity', 'distance_to_bone'])

    reg['SampleID'] = counter
    print(f"SampleID set to {counter} for {filename}")

        # Reorder columns
    new_column_order = ['x', 'y', 'label', 'cell_type', 'SampleID', 'cell_size']
    reg = reg[new_column_order]

        # Generate a new filename based on unique counter ID
    new_filename = f"patient{counter}_cell_positions.csv"

        # Save details to patient_file_data
    patient_file_data.append({
        "original_file_name": base_name,
        "patient_file_name": new_filename,
        "counter": counter
    })

        # Save all processed files to all_files_output_folder
    reg.to_csv(os.path.join(files_output_folder, new_filename), index=False)

        # Increment the counter for unique patient IDs
    counter += 1

print(f"Finished processing files")

# Save the patient file summary
summary_df = pd.DataFrame(patient_file_data)
summary_df.to_csv(os.path.join(files_output_folder, "patient_file_summary.csv"), index=False)

print("Patient file summary saved as patient_file_summary.csv")


# Load the patient_file_data CSV
patient_file_data = pd.read_csv("TMENS_analysis/data/cell_positions_data/patient_file_summary.csv")

# List to store data from all intensity files
all_data = []

# Directory containing the intensity files
intensities_directory = 'updated_csvs/intensities/' 

intensities_output_folder = "TMENS_analysis/data/cell_intensities_data/"
os.makedirs(intensities_output_folder, exist_ok=True)

# Iterate through each file in the intensities directory
for filename in os.listdir(intensities_directory):
    if filename.endswith(".csv"):
        # Construct the full file path
        file_path = os.path.join(intensities_directory, filename)
        
        # Check if the file exists in the patient_file_data DataFrame
        match = patient_file_data[patient_file_data['original_file_name'] == filename]
        
        if not match.empty:
            # Get the corresponding counter (SampleID) value
            sample_id = match.iloc[0]['counter']
            
            # Load the current intensities file
            intensities = pd.read_csv(file_path)
            # Add SampleID column with the matched counter value
            intensities['SampleID'] = sample_id

            
            intensities.rename(columns={'Object': 'cellLabelInImage'}, inplace=True)
        
            new_column_order = ['SampleID', 'cellLabelInImage', 'CD38', 'Perilipin', 'Vimentin', 'B4GALT1', 'MPO',
       'CathepsinK', 'ATP5A', 'RUNX2', 'HIF1A', 'CD11b', 'CD45', 'CS', 'CD11c',
       'CD36', 'CD4', 'CD34', 'CD68', 'IL32', 'IDO', 'CD8', 'GranzymeK',
       'PKM2', 'IRF4', 'GLUT1', 'GranzymeB', 'Ki67', 'CollagenTypeI', 'CD3',
       'CPT1A', 'CD98', 'HLA-DR', 'ST6GAL1', 'CD138']
        
            # Ensure only available columns are used
            intensities = intensities[[col for col in new_column_order if col in intensities.columns]]
                    # Save or process intensities_df as needed
                    # Optionally, overwrite the file or save to a new directory
                    # intensities_df.to_csv(file_path, index=False)
            # Add the DataFrame to all_data for merging later
            all_data.append(intensities)
            
            print(f"Processed file: {filename} with SampleID: {sample_id}")
        else:
            print(f"No match found for file: {filename}")


# Concatenate all DataFrames into one big DataFrame
merged_data = pd.concat(all_data, ignore_index=True)

# Save the merged DataFrame as 'cellData.csv'
merged_data.to_csv(os.path.join(intensities_output_folder, 'cellData.csv'), index=False)
print(f"All CSV files have been merged and saved to {os.path.join(intensities_output_folder, 'cellData.csv')}")
