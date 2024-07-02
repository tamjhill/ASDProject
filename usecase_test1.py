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

def process_csv(file_path, search_terms, search_column, return_column):
    results = []

    # Prepare search terms as whole word patterns
    search_patterns = [re.compile(r'\b' + re.escape(term) + r'\b', re.IGNORECASE) for term in search_terms]

    with open(file_path, 'r', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        
        for row in reader:
            search_text = row[search_column]
            return_value = row[return_column]

            # Check for matching terms
            matched_terms = set(term for term, pattern in zip(search_terms, search_patterns) if pattern.search(search_text))

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