import csv
import rdflib
from rdflib import Namespace, URIRef, Literal, BNode
import os
from rdflib.namespace import XSD
import requests

# Create the graph object which holds the triples
graph = rdflib.Graph()

# Define namespaces
BIOLINK = Namespace("https://w3id.org/biolink/vocab/")
ENSEMBL = Namespace("http://identifiers.org/ensembl/")
NCBIGENE = Namespace("http://identifiers.org/ncbigene/")
RDFS = Namespace("http://www.w3.org/2000/01/rdf-schema#")
RDF = Namespace("http://www.w3.org/1999/02/22-rdf-syntax-ns#")
SCHEMA = Namespace("https://schema.org/")
EDAM = Namespace("http://edamontology.org/")
DOI = Namespace("https://doi.org/")
DCT = Namespace("http://purl.org/dc/terms/")
PMC = Namespace("https://pubmed.ncbi.nlm.nih.gov/")
OWL = Namespace("http://www.w3.org/2002/07/owl#")


# Bind namespaces to prefixes
graph.bind("biolink", BIOLINK)
graph.bind("ensembl", ENSEMBL)
graph.bind("ncbigene", NCBIGENE)
graph.bind("rdfs", RDFS)
graph.bind("rdf", RDF)
graph.bind("schema", SCHEMA)
graph.bind("edam", EDAM)
graph.bind("doi", DOI)
graph.bind("dct", DCT)
graph.bind("pmc", PMC)
graph.bind("owl", OWL)

# Root directory to search for CSV files
root_dir = "C:\\Users\\tamjh\\CodeProjects\\ASDProject\\data"

def process_metadata_csv(csv_file_path):
    with open(csv_file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        
        for row in reader:
            if 'pmid' in row and row['pmid']:
                pmid_uri = PMC[row['pmid']]
                
                for column, value in row.items():
                    if column == 'doi' and value:
                        graph.add((pmid_uri, DCT.identifier, DOI[value]))
                for column, value in row.items():
                    if column == 'title' and value:
                        graph.add((pmid_uri, DCT.title, Literal(value)))
                for column, value in row.items():
                    if column == 'year' and value:
                        graph.add((pmid_uri, DCT.date, Literal(value)))
                for column, value in row.items():
                    if column == 'journal' and value:
                        graph.add((pmid_uri, DCT.publisher, Literal(value)))


def get_gene_id(gene_name):
    search_url = "https://api-v3.monarchinitiative.org/v3/api/search"
    params = {
        "q": gene_name,
        "category": "biolink:Gene",
        "limit": 1
    }
    try:
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        data = response.json()
        if data['items']:
            # Extract the HGNC ID from the 'id' field
            monarch_id = data['items'][0]['id']
            if monarch_id.startswith('HGNC:'):
                return f"https://monarchinitiative.org/{monarch_id}"
            # If it's not an HGNC ID, check if it's in the 'equivalent_identifiers'
            for eq_id in data['items'][0].get('equivalent_identifiers', []):
                if eq_id.startswith('HGNC:'):
                    return f"https://monarchinitiative.org/{eq_id}"
    except requests.exceptions.RequestException as e:
        print(f"An error occurred while searching for {gene_name}: {e}")
    return None

def process_regular_csv(csv_file_path):
    filename = os.path.splitext(os.path.basename(csv_file_path))[0]
    pmid = os.path.basename(os.path.dirname(csv_file_path))
    
    dataset_uri = BNode()
    pmid_uri = PMC[pmid]
    
    graph.add((dataset_uri, RDF.type, BIOLINK.Dataset))
    graph.add((pmid_uri, EDAM.has_output, dataset_uri))
    graph.add((dataset_uri, OWL.sameAs, Literal(filename)))
    graph.add((pmid_uri, RDF.type, DCT.identifier))
    
    possible_gene_names = ['symbol', 'genesymbol', 'genename', 'geneid', 'ensembl', 'ensemblid', 'entrez', 'entrezid', 'ncbi', 'ncbiid']
    possible_log_names = ['log2', 'lf2', 'lfc2', 'logfold2', 'log2fc', 'logfoldchange', 'logfold', 'lf', 'logfc', 'foldchange', 'fc', 'lfc', 'enrichment']
    possible_pval_names = ['padj', 'adjp', 'pvalueadj', 'adjpvalue', 'pvaladj', 'adjpval', 'pvadj', 'adjpv', 'pvalue', 'pval', 'pv']
    
    with open(csv_file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            row_node = BNode()
            graph.add((dataset_uri, EDAM.has_output, row_node))
            graph.add((row_node, RDF.type, EDAM.data))
            
            gene_added = log_added = pval_added = False
            
            for column, value in row.items():
                if value:  # Only process non-empty values
                    column_lower = column.lower()
                    
                    # Check for gene names
                    if any(name in column_lower for name in possible_gene_names):
                        if not gene_added:
                            monarch_uri = get_gene_id(value)
                            if monarch_uri:
                                graph.add((row_node, BIOLINK.Gene, URIRef(monarch_uri)))
                                gene_added = True
                            else:
                                if 'ensembl' in column_lower:
                                    graph.add((row_node, ENSEMBL.id, Literal(value)))
                                elif 'entrez' in column_lower or 'ncbi' in column_lower:
                                    graph.add((row_node, NCBIGENE.id, Literal(value)))
                                else:
                                    graph.add((row_node, BIOLINK.symbol, Literal(value)))
                                gene_added = True
                    
                    # Check for log fold change
                    elif any(name in column_lower for name in possible_log_names):
                        if not log_added:
                            graph.add((row_node, EDAM.data_3754, Literal(value)))
                            log_added = True
                        else:
                            predicate = URIRef(f"rdf:predicate/{column}")
                            graph.add((row_node, predicate, Literal(value)))
                    
                    # Check for p-values
                    elif any(name in column_lower for name in possible_pval_names):
                        if not pval_added:
                            graph.add((row_node, EDAM.data_2082, Literal(value)))
                            pval_added = True
                        else:
                            predicate = URIRef(f"rdf:predicate/{column}")
                            graph.add((row_node, predicate, Literal(value)))
                    
                    # Add any other columns as generic predicates
                    else:
                        predicate = URIRef(f"rdf:predicate/{column}")
                        graph.add((row_node, predicate, Literal(value)))


for dirpath, dirnames, filenames in os.walk(root_dir):
    for filename in filenames:
        if filename.endswith('.csv'):
            csv_file_path = os.path.join(dirpath, filename)
            print(f"Processing file: {csv_file_path}")
            
            if filename == 'asd_article_metadata.csv':
                process_metadata_csv(csv_file_path)
            else:
                process_regular_csv(csv_file_path)

graph.serialize(destination='main_graph.nt', format='nt', encoding= "utf-8" )
print("Combined graph has been serialized to main_graph.nt")
