import pandas as pd
import os

root = "data"
filename1 = "asd_article_metadata.csv"
filename2 = "dataset_test1.csv"

a = pd.read_csv(os.path.join(root, filename1))
b = pd.read_csv(os.path.join(root, filename2))
#a = b.dropna(axis=1)
merged = a.merge(b, on='pmid')
merged.to_csv("output.csv", index=False)