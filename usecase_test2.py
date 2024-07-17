import spacy
import csv
import re

def get_csv_column(file_path, column_name):
    data = []
    with open(file_path, 'r') as csvfile:
        csvreader = csv.DictReader(csvfile)
        for row in csvreader:
            data.append(row[column_name])
    return data

#print(len(column_data))

import csv
import spacy

def process_csv(file_path, search_terms, search_column, return_column):
    nlp = spacy.load("en_core_web_sm")

    results = []
    search_terms_lower = [term.lower() for term in search_terms]

    with open(file_path, 'r', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        
        for row in reader:
            search_text = row[search_column]
            return_value = row[return_column]

            doc = nlp(search_text)

            matched_terms = set(token.text for token in doc 
                                if token.text.lower() in search_terms_lower 
                                and not token.is_punct
                                and token.pos_ != "VERB")

            if matched_terms:
                results.append({
                    'matched_terms': list(matched_terms),
                    'return_value': return_value
                })

    return results

abstract_file = 'data/asd_article_metadata.csv'

gene_file = 'gene_list.csv'
column_name = 'Gene name'
gene_list = get_csv_column(gene_file, column_name)
search_column = 'abstract'
return_column = 'doi'

results = process_csv(abstract_file, gene_list, search_column, return_column)

for result in results:
    print(f"Matched terms: {result['matched_terms']}")
    print(f"Return value: {result['return_value']}")
    print()

"""returns:
Matched terms: ['set']
Return value: 10.3389/fncir.2022.982721

Matched terms: ['CD14']
Return value: 10.33696/immunology.4.146

Matched terms: ['impact', 'TEs', 'ATRX', 'set']
Return value: 10.1186/s13229-023-00554-5

Matched terms: ['GRN']
Return value: 10.1126/science.adh2602
"""