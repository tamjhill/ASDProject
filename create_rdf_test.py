import csv
import rdflib
from rdflib import Namespace, URIRef, Literal, BNode
import os
from rdflib.namespace import XSD

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
PMC = Namespace("https://www.ncbi.nlm.nih.gov/pmc/articles/PMID")
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

def process_regular_csv(csv_file_path):
    # Get the filename without extension and the folder name (which will be used as the DOI)
    filename = os.path.splitext(os.path.basename(csv_file_path))[0]
    pmid = os.path.basename(os.path.dirname(csv_file_path))
    #doi = folder_name.replace("_", "/")
    
    # Create a URI for the dataset and for the DOI
    dataset_uri = BNode()
    pmid_uri = PMC[pmid]
    
    # Add dataset type and the triple linking DOI to dataset
    graph.add((dataset_uri, RDF.type, BIOLINK.Dataset))
    graph.add((pmid_uri, EDAM.has_output, dataset_uri))
    graph.add((dataset_uri, OWL.sameAs, Literal(filename)))
    graph.add((pmid_uri, RDF.type, DCT.identifier))
    
    with open(csv_file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        
        for row in reader:
            # Create a blank node for each row
            row_node = BNode()
            graph.add((dataset_uri, EDAM.has_output, row_node))
            graph.add((row_node, RDF.type, EDAM.data))
            possible_column_names = ['symbol', 'gene', 'gene_symbol', 'gene_name', 'genesymbol', 'genename']
    
            for column, value in row.items():
                if value:  # Only add triples if the value is not empty
                    if 'logfold' in column.lower():
                        graph.add((row_node, EDAM.data_3754, Literal(value)))
                    elif 'pvalue' in column.lower():
                        graph.add((row_node, EDAM.data_2082, Literal(value)))
                    elif column.lower() in possible_column_names :
                        graph.add((row_node, BIOLINK.symbol, Literal(value)))
                    elif 'ensembl' in column.lower():
                        graph.add((row_node, ENSEMBL.id, Literal(value)))
                    elif ('entrez' in column.lower()) or ('ncbi' in column.lower()):
                        graph.add((row_node, NCBIGENE.id, Literal(value)))
                    else:
                        continue
#                    else:
#                        # For all other columns, use a generic predicate
#                        predicate = URIRef(f"rdf:predicate/{column}")
#                        graph.add((row_node, predicate, Literal(value)))

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
