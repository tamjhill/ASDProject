from eparse.core import get_df_from_file
import os
import pandas as pd

folderpath = "data/supp_data/10.1038_s41586-022-05377-7"
csv = 'csv'

# Get the list of tables
table_list = list(get_df_from_file('data/supp_data/10.1038_s41586-022-05377-7/41586_2022_5377_MOESM8_ESM.xlsx'))

# Iterate through the tables
for i, table in enumerate(table_list, 1):
    # Check if the table is a tuple and extract the DataFrame
    if isinstance(table, tuple):
        df = table[0]  # Assuming the DataFrame is the first element of the tuple
    elif isinstance(table, pd.DataFrame):
        df = table
    else:
        print(f"Unexpected type for table {i}: {type(table)}")
        continue

    # Save the DataFrame to CSV
    filename = f"df_{i}.{csv}"
    filepath = os.path.join(folderpath, filename)
    df.to_csv(filepath, index=False)