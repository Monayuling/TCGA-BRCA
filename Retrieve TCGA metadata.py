#Retrieve Clinical Data including ER, PR, HER2 from GDC Portal


# Install and load TCGAbiolinks package
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)

clinical_query <- GDCquery(  
	project = "TCGA-BRCA",   
	data.category = "Clinical",  
	data.type = "Clinical Supplement",   
	data.format = "BCR XML"
)

GDCdownload(clinical_query)
clinical_tab_all <- GDCprepare(clinical_query)


#New code based on one xml file
#Read all xml files and from GDC Portal and compile into large matrix (1000+ xml files with one sample in formation in one file)
import os
import xml.etree.ElementTree as ET
import pandas as pd

# Define the path to the root directory containing all the XML files
root_directory = '/content/drive/MyDrive/Clinical_Supplement'

# Create an empty list to store DataFrames for each XML file
dataframes = []

# List of desired text elements
desired_text = [
    "bcr_patient_barcode", 
    "gender",
    "tissue_source_site", 
    "history_of_neoadjuvant_treatment", 
    "year_of_initial_pathologic_diagnosis", 
    "menopause_status", 
    "estrogen_receptor_status", 
    "progesterone_receptor_status", 
    "her2_new_chromosone_17_signal_ratio_value", 
    "her2_immunohistochemistry_level_result", 
    "er_level_cell_percentage_category", 
    "progesterone_receptor_level_cell_percent_category", 
    "her2_new_immunohistochemistry_receptor_status", 
    "distant_metastasis_present_ind2", 
    "pathologic_stage", 
    "pathologic_T", 
    "pathologic_N", 
    "pathologic_M", 
    "new_tumor_event_after_initial_treatment", 
    "radiation_therapy", 
    "vital_status", 
    "days_to_last_followup", 
    "days_to_death", 
    "person_neoplasm_cancer_status", 
    "days_to_new_tumor_event_after_initial_treatment", 
    "new_neoplasm_event_occurrence_anatomic_site", 
    "year_of_form_completion", 
    "therapy_type", 
    "regimen_indication", 
    "age_at_initial_pathologic_diagnosis", 
    "year_of_initial_pathologic_diagnosis"
]

# Function to extract data from an XML file and return a DataFrame
def process_xml_file(xml_file_path):
    data = {}
    tree = ET.parse(xml_file_path)
    root = tree.getroot()
    for element in root.iter():
        tag = element.tag.split('}')[1]
        data[tag] = element.text
    df = pd.DataFrame([data])
    return df

# Iterate through the root directory and its subdirectories
for directory_path, subdirectories, files in os.walk(root_directory):
    for file_name in files:
        if file_name.endswith('.xml'):
            file_path = os.path.join(directory_path, file_name)
            # Process the XML file and append the DataFrame to the list
            df = process_xml_file(file_path)
            dataframes.append(df)

# Concatenate all DataFrames into a single DataFrame
final_df = pd.concat(dataframes, ignore_index=True)

# Filter columns based on desired text elements
filtered_columns = [col for col in final_df.columns if any(text in col for text in desired_text)]
new_df = final_df[filtered_columns]

# Display the resulting DataFrame
new_df

#------------------------------------------------
# Filter out male patients and patients who had undergone neoadjuvant chemotherapy
filtered_df = new_df[(new_df['gender'] != 'MALE') & (new_df['history_of_neoadjuvant_treatment'] == 'No')]

# Reset the index of the filtered DataFrame if needed
filtered_df.reset_index(drop=True, inplace=True)

# Display the filtered DataFrame
filtered_df

# Specify the file path for the output text file
filtered_file_path = '/content/drive/MyDrive/filtered.txt' #Save file as 'filtered.txt'

# Extract the 'patient_bar_code' column and write it to a text file
filtered_df['bcr_patient_barcode'].to_csv(filtered_file_path, header=False, index=False)

#------------------------------------------------------
#Get desired subgroups, and create a txt file with patient bar codes of selected patients
# Get parameter names 
new_df.columns

#Get more concise dataframe
required = ['bcr_patient_barcode', 'tissue_source_site','age_at_initial_pathologic_diagnosis','year_of_initial_pathologic_diagnosis', 'person_neoplasm_cancer_status',
       'menopause_status', 'breast_carcinoma_progesterone_receptor_status',
       'her2_immunohistochemistry_level_result',
       'breast_carcinoma_estrogen_receptor_status',
       'er_level_cell_percentage_category',
       'progesterone_receptor_level_cell_percent_category',
       'distant_metastasis_present_ind2','pathologic_stage', 'pathologic_T', 'pathologic_N', 'pathologic_M','days_to_new_tumor_event_after_initial_treatment',
       'new_neoplasm_event_occurrence_anatomic_site']

c_columns = [col for col in filtered_df.columns if any(text in col for text in required)]
c_df = filtered_df[c_columns]
c_df

#Get luminal subtype
luminal = filtered_df[(filtered_df['breast_carcinoma_estrogen_receptor_status'] == 'Positive') | (filtered_df['breast_carcinoma_progesterone_receptor_status'] == 'Positive')]
luminal

# Specify the file path for the output text file
luminal_file_path = '/content/drive/MyDrive/luminal.txt'

# Extract the 'patient_bar_code' column and write it to a text file
luminal['bcr_patient_barcode'].to_csv(luminal_file_path, header=False, index=False)

#Get pre and post menopausal patients
#Premenopausal
premenopausal = luminal[(filtered_df['menopause_status'] != 'Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)')]
premenopausal

# Specify the file path for the output text file
premenopausal_file_path = '/content/drive/MyDrive/premenopausal.txt'

# Extract the 'patient_bar_code' column and write it to a text file
premenopausal['bcr_patient_barcode'].to_csv(premenopausal_file_path, header=False, index=False)

#Postmenopausal
postmenopausal = luminal[(filtered_df['menopause_status'] == 'Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)')]
postmenopausal

# Specify the file path for the output text file
postmenopausal_file_path = '/content/drive/MyDrive/postmenopausal.txt'

# Extract the 'patient_bar_code' column and write it to a text file
postmenopausal['bcr_patient_barcode'].to_csv(postmenopausal_file_path, header=False, index=False)

#Find young and old patients
#Convert the 'age_at_initial_pathologic_diagnosis' column to integers
luminal['age_at_initial_pathologic_diagnosis'] = pd.to_numeric(luminal['age_at_initial_pathologic_diagnosis'], errors='coerce')

# Filter the 'luminal' DataFrame for individuals with age less than or equal to 40
very_young = luminal[luminal['age_at_initial_pathologic_diagnosis'] <= 40]
very_young_file_path = '/content/drive/MyDrive/very_young.txt'
very_young['bcr_patient_barcode'].to_csv(very_young_file_path, header=False, index=False)

# Filter the 'luminal' DataFrame for individuals with age older than or equal to 80
very_old = luminal[luminal['age_at_initial_pathologic_diagnosis'] >= 80]
very_old_file_path = '/content/drive/MyDrive/very_old.txt'
very_old['bcr_patient_barcode'].to_csv(very_old_file_path, header=False, index=False)

# Extract young patients: <= 50 years old
# Filter the 'luminal' DataFrame for individuals with age less than or equal to 50
young = luminal[luminal['age_at_initial_pathologic_diagnosis'] <= 50]
young_file_path = '/content/drive/MyDrive/young.txt'
young['bcr_patient_barcode'].to_csv(young_file_path, header=False, index=False)

# Filter the 'luminal' DataFrame for individuals with age older than 50
old = luminal[luminal['age_at_initial_pathologic_diagnosis'] > 50]
old_file_path = '/content/drive/MyDrive/old.txt'
old['bcr_patient_barcode'].to_csv(old_file_path, header=False, index=False)

#---------------------------------------
#Create expression matrices with desired parameters:

import pandas as pd

# Define the file paths
excel_file_path = '/content/drive/MyDrive/exp_matrix.xlsx'  # Path to the Excel file

def extract(parameter):
  # Read the patient bar codes from the text file into a list
  text_file_path = '/content/drive/MyDrive/' + parameter + '.txt'  # Path to the text file with patient bar codes
  with open(text_file_path, 'r') as file:
    patient_barcodes = [line.strip() for line in file]

  # Read the expression matrix from the Excel file
  expression_matrix = pd.read_excel(excel_file_path)

  # Use regular expressions to match column names based on the first portion of patient barcodes
  matching_columns = expression_matrix.filter(regex=f'^({"|".join(patient_barcodes)}).*', axis=1)

  # Save the extracted columns to a new text file
  output_text_file_path = '/content/drive/MyDrive/' + parameter + '_expression.txt'
  matching_columns.to_csv(output_text_file_path, sep='\t', index=False)

  print("Columns extracted and saved to", output_text_file_path)

extract('luminal')
extract('very_young')
extract('very_old')
extract('premenopausal')
extract('postmenopausal')
extract('young')
extract('old')
