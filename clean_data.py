#import excel
#clean data


import pandas as pd
import os

root = 'data'
first_dir = 'supp_data'
second_dir = 's00335-024-10036-5'
filename = '335_2024_10036_MOESM2_ESM.xlsx'
file_path = os.path.join(root, first_dir, second_dir, filename)

df = pd.read_excel(file_path)
print(df)





# Iterate over files in a directory (inc. sub-directories)
#data_dir = 'data'
#for root, dirs, files in os.walk(data_dir):
#    for filename in files:
#        file_path = os.path.join(root, filename)
#        
#        if filename.endswith('.xlsx', '.xls'):
#            df = pd.read_excel()