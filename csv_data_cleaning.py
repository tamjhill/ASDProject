import pandas as pd
import os

def process_csv_file(file_path):
    # Read the CSV file
    df = pd.read_csv(file_path)
    
    # Convert column titles to lowercase and replace spaces with underscores
    df.columns = [col.lower().replace(' ', '_') for col in df.columns]
    
    # Remove rows with multiple NAs or blanks
<<<<<<< HEAD
=======
    # This will keep rows where at least half of the columns have non-null values
>>>>>>> 85d51d638d226caa742fc1b6f2febc20341adfc5
    df_cleaned = df.dropna(thresh=len(df.columns)//5)
    
    # Save the processed dataframe back to CSV
    df_cleaned.to_csv(file_path, index=False)
    print(f"Processed and overwritten: {file_path}")

def process_data_folder(data_folder):
<<<<<<< HEAD
=======
    # Walk through all subdirectories
>>>>>>> 85d51d638d226caa742fc1b6f2febc20341adfc5
    for root, dirs, files in os.walk(data_folder):
        for file in files:
            if file.endswith('.csv'):
                file_path = os.path.join(root, file)
                print(f"Processing file: {file_path}")
                process_csv_file(file_path)

<<<<<<< HEAD
=======
# Usage
>>>>>>> 85d51d638d226caa742fc1b6f2febc20341adfc5
data_folder = "data\supp_data"
process_data_folder(data_folder)