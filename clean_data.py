#import excel
#clean data


import pandas as pd
import os

#root = 'data'
#first_dir = 'supp_data'
#second_dir = 's00335-024-10036-5'
#filename = '335_2024_10036_MOESM2_ESM.xlsx'
#file_path = os.path.join(root, first_dir, second_dir, filename)
#
#df = pd.read_excel(file_path)
#print(df)


import csv

# Set the input and output file paths
input_file = 'data\kg.csv'
output_file = 'data\kg_new.csv'

# Define the list of relation keywords to filter out
relation_keywords = ['contraindication', 'drug_protein', 'drug_drug', 'drug_effect', 
            'cellcomp_cellcomp', 'exposure_protein', 'exposure_exposure', 
            'exposure_molfunc', 'exposure_cellcomp']

# Define the list of keywords to filter out in 'x_source' column
source_keywords = ['DrugBank', 'CTD']

# Read the input file
with open(input_file, 'r', encoding='utf-8') as file:
    reader = csv.DictReader(file)
    data = list(reader)

filtered_data = [row for row in data
                 if row['relation'].lower() not in [keyword.lower() for keyword in relation_keywords]
                 and row['x_source'].lower() not in [keyword.lower() for keyword in source_keywords]]

# Write the filtered data to the output file
with open(output_file, 'w', newline='', encoding='utf-8') as file:
    fieldnames = data[0].keys()
    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(filtered_data)

print(f"Filtered data saved to {output_file}")


# Iterate over files in a directory (inc. sub-directories)
#data_dir = 'data'
#for root, dirs, files in os.walk(data_dir):
#    for filename in files:
#        file_path = os.path.join(root, filename)
#        
#        if filename.endswith('.xlsx', '.xls'):
#            df = pd.read_excel()


#import csv
#def count_rows(file_path):
#    with open(file_path, 'r', encoding='utf-8') as file:
#        reader = csv.reader(file)
#        row_count = sum(1 for row in reader)
#        return row_count


#input_file = 'data\edges.csv'
#output_file = 'data\edges_new.csv'
#num_rows = count_rows(input_file)
#print(f"The CSV file '{input_file}' has {num_rows} rows.")

#num_rows_new = count_rows(output_file)
#print(f"The CSV file '{output_file}' has {num_rows_new} rows.")

import csv

# Set the input and output file paths
input_file = 'data\kg_GO.csv'
output_file = 'data\kg_GO_genes.csv'

# Define the keyword to filter for in the 'y_source' column
keyword = 'gene/protein'

# Read the input file
with open(input_file, 'r', encoding='utf-8') as file:
    reader = csv.DictReader(file)
    data = list(reader)

# Filter out rows where the 'y_source' column does not contain the keyword
filtered_data = [row for row in data if keyword.lower() in row['x_type'].lower()]

# Write the filtered data to the output file
with open(output_file, 'w', newline='', encoding='utf-8') as file:
    fieldnames = data[0].keys()
    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(filtered_data)

print(f"Filtered data saved to {output_file}")