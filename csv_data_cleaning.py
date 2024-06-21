import pandas as pd
import os

def process_csv_file(file_path):
    # Read the CSV file
    df = pd.read_csv(file_path)
    
    # Convert column titles to lowercase and replace spaces with underscores
    df.columns = [col.lower().replace(' ', '_') for col in df.columns]
    
    # Remove rows with multiple NAs or blanks
    # This will keep rows where at least half of the columns have non-null values
    df_cleaned = df.dropna(thresh=len(df.columns)//5)
    
    # Save the processed dataframe back to CSV
    df_cleaned.to_csv(file_path, index=False)
    print(f"Processed and overwritten: {file_path}")

def process_data_folder(data_folder):
    # Walk through all subdirectories
    for root, dirs, files in os.walk(data_folder):
        for file in files:
            if file.endswith('.csv'):
                file_path = os.path.join(root, file)
                print(f"Processing file: {file_path}")
                process_csv_file(file_path)

# Usage
data_folder = "data\supp_data"
process_data_folder(data_folder)