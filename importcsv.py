import csv

with open('data/asd_article_metadata.csv', 'r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        print(row['pmid'])