"""takes excel files retrieved from pubmed, checks if they contain a column re. log fold changes, and saves sheet as csv file"""

import pandas as pd
import os
from openpyxl import load_workbook
import xlrd
import csv

def process_excel_file(file_path):
    # Load the Excel file
    wb = load_workbook(filename=file_path, read_only=True)
    
    # Get the directory of the input file
    output_dir = os.path.dirname(file_path)
    
    # Iterate through each sheet
    for sheet_name in wb.sheetnames:
        df = pd.read_excel(file_path, sheet_name=sheet_name)
        
        # Check if the sheet has a column with "log fold change" or similar
        log_fold_col = None
        for col in df.columns:
            if any(phrase in col.lower().replace('-', ' ') for phrase in ['log fold change', 'log fold', 'log fold2', 'lf', 'enrichment', 'logfc', 'fold change', 'fc', 'log2']):
                log_fold_col = col
                break
        
        # If a matching column is found, save the sheet as CSV
        if log_fold_col:
            newfile = sheet_name.replace(" ", "")
            output_file = os.path.join(output_dir, f"{newfile}.csv")
            df.to_csv(output_file, index=False)
            print(f"Saved {sheet_name} as CSV: {output_file}")
        else:
            print(f"Skipped {sheet_name} in {file_path}: No 'log fold change' column found")

def process_old_file(file_path):
    # Load the xls file
    wb = xlrd.open_workbook(file_path)

    output_dir = os.path.dirname(file_path)
    
    # Iterate through each sheet
    for sheet in wb.sheets():
        print(f"Processing sheet: {sheet.name}")

        # Get the headers
        headers = [sheet.cell_value(0, col) for col in range(sheet.ncols)]

        # Check if the sheet has a column with "log fold change" or similar
        log_fold_col = None
        for col, header in enumerate(headers):
            if any(phrase in header.lower().replace('-', ' ') for phrase in ['log fold change', 'log fold', 'log fold2', 'lf', 'enrichment', 'logfc', 'fold change', 'fc', 'log2']):
                log_fold_col = col
                break

        # If a matching column is found, save the sheet as CSV
        if log_fold_col is not None:
            newfile = sheet.name.replace(" ", "")
            output_file = os.path.join(output_dir, f"{newfile}.csv")

            with open(output_file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter=",")
                writer.writerow(headers)  # Write headers

                for row in range(1, sheet.nrows):  # Start from 1 to skip header
                    row_data = [sheet.cell_value(row, col) for col in range(sheet.ncols)]
                    writer.writerow(row_data)

            print(f"Saved {sheet.name} as CSV: {output_file}")
        else:
            print(f"Skipped {sheet.name} in {file_path}: No 'log fold change' column found")

def process_data_folder(data_folder):
    # Walk through all subdirectories
    for root, dirs, files in os.walk(data_folder):
        for file in files:
            if file.endswith('.xlsx'):
                file_path = os.path.join(root, file)
                print(f"Processing file: {file_path}")
                process_excel_file(file_path)
            elif file.endswith('.xls'):
                file_path = os.path.join(root, file)
                print(f"Processing file: {file_path}")
                process_old_file(file_path)

data_folder = "data\\supp_data"
process_data_folder(data_folder)