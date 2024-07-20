import pandas as pd
import os
import re

def process_csv_file(file_path):
    # Read the CSV file
    df = pd.read_csv(file_path)
    # Remove rows with multiple NAs or blanks
    df_cleaned = df.dropna(thresh=len(df.columns)//4)
    # Process column names
    df_cleaned.columns = df_cleaned.columns.map(lambda x: re.sub(r'[\s\-_]', '', x.lower()))
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