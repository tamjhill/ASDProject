import pandas as pd
import os
import re

def process_csv_file(file_path):
    # Read the CSV file
    df = pd.read_csv(file_path)
    
    # Function to process column names
    def process_column_name(col):
        col = col.lower().replace(' ', '_')
        if re.search(r'p[-_]?val(ue)?', col):
            return 'pvalue'
        elif re.search(r'(log[-_]?fold[-_]?(change|2)?|lf2?|fold[-_]?change)', col):
            return 'logfold'
        return col
    
    # Apply the processing function to column names
    df.columns = [process_column_name(col) for col in df.columns]
    # Remove rows with multiple NAs or blanks
    df_cleaned = df.dropna(thresh=len(df.columns)//4)
    # Save the processed dataframe back to CSV
    df_cleaned.to_csv(file_path, index=False, sep=",")
    print(f"Processed and overwritten: {file_path}")

def process_data_folder(data_folder):
    for root, dirs, files in os.walk(data_folder):
        for file in files:
            if file.endswith('.csv'):
                file_path = os.path.join(root, file)
                print(f"Processing file: {file_path}")
                process_csv_file(file_path)

data_folder = "data\supp_data"
process_data_folder(data_folder)