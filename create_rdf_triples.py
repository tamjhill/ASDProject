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
MONARCH = Namespace("https://monarchinitiative.org/")

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
graph.bind("monarch", MONARCH)


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
    genefile = 'gene_ids.txt'
    #print(f"Searching for gene: {gene_name}")
    with open(genefile, 'r') as file:
        for line in file:
            if gene_name.upper() in line.upper():
                #print(f"{gene_name} HGNC ID found")
                return line.split('\t')[0]
    #print(f"HGNC ID not found for {gene_name}")
    return None
    
    
def process_regular_csv(csv_file_path, matched_genes, unmatched_genes):
    filename = os.path.splitext(os.path.basename(csv_file_path))[0]
    pmid = os.path.basename(os.path.dirname(csv_file_path))
    
    dataset_uri = BNode()
    pmid_uri = PMC[pmid]
    
    graph.add((dataset_uri, RDF.type, BIOLINK.Dataset))
    graph.add((pmid_uri, EDAM.has_output, dataset_uri))
    graph.add((dataset_uri, OWL.sameAs, Literal(filename)))
    graph.add((pmid_uri, RDF.type, DCT.identifier))
    
    possible_gene_names = ['ensembl', 'symbol', 'genesymbol', 'genename', 'geneid', 'entrez', 'ncbi']
    possible_log_names = ['log2', 'lf2', 'lfc2', 'logfold2', 'log2fc', 'logfoldchange', 'logfold', 'lf', 'logfc', 'foldchange', 'fc', 'lfc', 'enrichment']
    possible_pval_names = ['padj', 'adjp', 'pvalueadj', 'adjpvalue', 'pvaladj', 'adjpval', 'pvadj', 'adjpv', 'pvalue', 'pval', 'pv']
    
    matched_genes = 0
    unmatched_genes = 0

    with open(csv_file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            row_node = BNode()
            graph.add((dataset_uri, EDAM.has_output, row_node))
            graph.add((row_node, RDF.type, EDAM.data))
            
            gene_added = log_added = pval_added = False
            
            # Try to find and process gene information
            for gene_column in possible_gene_names:
                for column, value in row.items():
                    if gene_column in column.lower() and value and not gene_added:
                        #print("gene column match")
                        monarch_uri = get_gene_id(value)
                        if monarch_uri:
                            graph.add((row_node, BIOLINK.Gene, MONARCH[monarch_uri]))
                            gene_added = True
                            matched_genes[0] += 1
                            #print(f"Added gene: {value} with HGNC ID: {monarch_uri}")
                            break
                        else:
                            if 'ensembl' in column.lower():
                                graph.add((row_node, ENSEMBL.id, Literal(value)))
                            elif 'entrez' in column.lower() or 'ncbi' in column.lower():
                                graph.add((row_node, NCBIGENE.id, Literal(value)))
                            else:
                                graph.add((row_node, BIOLINK.symbol, Literal(value)))
                            gene_added = True
                            unmatched_genes[0] += 1
                            #print(f"Added gene {value} without HGNC ID")
                            break
                if gene_added:
                    break
            
            # Then process the rest of the columns
            for column, value in row.items():
                if value:  # Only process non-empty values
                    column_lower = column.lower()
                    
                    # Check for log fold change
                    if any(name in column_lower for name in possible_log_names) and not log_added:
                        graph.add((row_node, EDAM.data_3754, Literal(value)))
                        log_added = True
                    
                    # Check for p-values
                    elif any(name in column_lower for name in possible_pval_names) and not pval_added:
                        graph.add((row_node, EDAM.data_2082, Literal(value)))
                        pval_added = True
                    
                    # Add any other columns as generic predicates
                    else:
                        predicate = URIRef(f"rdf:predicate/{column}")
                        graph.add((row_node, predicate, Literal(value)))
    return matched_genes, unmatched_genes

# Root directory to search for CSV files
root_dir = "C:\\Users\\tamjh\\CodeProjects\\ASDProject\\test_data"
matched_genes = [0]  
unmatched_genes = [0]  
total_files_processed = 0

for dirpath, dirnames, filenames in os.walk(root_dir):
    for filename in filenames:
        if filename.endswith('.csv'):
            csv_file_path = os.path.join(dirpath, filename)
            print(f"Processing file: {csv_file_path}")
            
            if filename == 'asd_article_metadata.csv':
                process_metadata_csv(csv_file_path)
            else:
                matched_genes, unmatched_genes = process_regular_csv(csv_file_path, matched_genes, unmatched_genes)
            total_files_processed += 1

print("\nProcessing complete. Summary:")
print(f"Total CSV files processed: {total_files_processed}")
print(f"Total matched genes: {matched_genes[0]}")
print(f"Total unmatched genes: {unmatched_genes[0]}")
print(f"Total genes processed: {matched_genes[0] + unmatched_genes[0]}")
graph.serialize(destination='main_graph.nt', format='nt', encoding= "utf-8" )
print("Combined graph has been serialized to test_graph.nt")


