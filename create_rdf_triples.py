import requests
import os
import csv
from rdflib import Graph, BNode, URIRef, Literal, Namespace
from functools import lru_cache

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


@lru_cache(maxsize=600)
def get_gene_id(gene_identifier, id_type=None):
    search_url = "https://api-v3.monarchinitiative.org/v3/api/search"
    
    # Determine the search query based on the ID type
    if id_type == 'ensembl':
        query = f"ENSEMBL:{gene_identifier}"
    elif id_type == 'entrez':
        query = f"NCBIGene:{gene_identifier}"
    else:
        query = gene_identifier

# set query parameters - for genes belonging to the human taxon
    params = {
        "q": query,
        "category": "biolink:Gene",
        "limit": 1,
        "in_taxon_label": "Homo sapiens"
    }
    try:
        #print(f"Searching for: {query}")
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        data = response.json()
        #print(f"API Response: {data}")
        if data['items']:
            for item in data['items']:
                #print(f"Checking item: {item}")
                if 'HGNC:' in item['id']:
                    #print("HGNC ID found in main ID")
                    return f"https://monarchinitiative.org/{item['id']}"
            # If no HGNC ID found in the main ID, check equivalent identifiers
            for eq_id in data['items'][0].get('equivalent_identifiers', []):
                #print(f"Checking equivalent ID: {eq_id}")
                if eq_id.startswith('HGNC:'):
                    #print("HGNC ID found in equivalent identifiers")
                    return f"https://monarchinitiative.org/{eq_id}"
        #else:
            #print("No items found in the API response")
    except requests.exceptions.RequestException as e:
        print(f"An error occurred while searching for {query}: {e}")
    print("No M-ID retrieved")
    return None

def get_gene_id_prep(ensembl_id, gene_symbol):
    result = get_gene_id(gene_symbol.upper())
    if result is None:
        result = get_gene_id(ensembl_id.upper(), id_type='ensembl')
    return result

# Test with MAD1L1
#result = get_gene_id_robust("ENSG00000002822", "MAD1L1")
#print(f"MAD1L1 Monarch ID: {result}")
#result = get_gene_id("MAD1L1")
#print(f"Final result: {result}")

def setup_graph(csv_file_path, graph):
    filename = os.path.splitext(os.path.basename(csv_file_path))[0]
    pmid = os.path.basename(os.path.dirname(csv_file_path))
    
    dataset_uri = BNode()
    pmid_uri = PMC[pmid]
    
    graph.add((dataset_uri, RDF.type, BIOLINK.Dataset))
    graph.add((pmid_uri, EDAM.has_output, dataset_uri))
    graph.add((dataset_uri, OWL.sameAs, Literal(filename)))
    graph.add((pmid_uri, RDF.type, DCT.identifier))
    #print("Graph initialised")
    return dataset_uri

def process_gene(row_node, column, value, gene_added, ensembl_id):
    if not gene_added:
        ensembl_id = next((row[col] for col in row if any(name in col.lower() for name in ['ensembl', 'geneid', 'ensemblid']) and row[col]), '')
        monarch_uri = get_gene_id_prep(ensembl_id, value)
        if monarch_uri:
            return [(row_node, BIOLINK.Gene, URIRef(monarch_uri))], True
        elif 'ensembl' in column.lower():
            return [(row_node, ENSEMBL.id, Literal(value))], True
        elif 'entrez' in column.lower() or 'ncbi' in column.lower():
            return [(row_node, NCBIGENE.id, Literal(value))], True
        else:
            return [(row_node, BIOLINK.symbol, Literal(value))], True
    return [], gene_added

def process_log_fold(row_node, column, value, log_added):
    if not log_added:
        #print("log-fold data added to the graph")
        return [(row_node, EDAM.data_3754, Literal(value))], True
    else:
        predicate = URIRef(f"rdf:predicate/{column}")
        return [(row_node, predicate, Literal(value))], log_added

def process_pvalue(row_node, column, value, pval_added):
    if not pval_added:
        #print("pvalue data added to the graph")
        return [(row_node, EDAM.data_2082, Literal(value))], True
    else:
        predicate = URIRef(f"rdf:predicate/{column}")
        return [(row_node, predicate, Literal(value))], pval_added

def process_other(row_node, column, value):
    predicate = URIRef(f"rdf:predicate/{column}")
    return [(row_node, predicate, Literal(value))]

def process_row(row, dataset_uri):
    row_node = BNode()
    triples = [(dataset_uri, EDAM.has_output, row_node),
               (row_node, RDF.type, EDAM.data)]
    
    gene_added = log_added = pval_added = False

    possible_gene_names = ['symbol', 'genesymbol', 'genename', 'geneid', 'ensembl', 'ensemblid', 'entrez', 'entrezid', 'ncbi', 'ncbiid']
    possible_log_names = ['log2', 'lf2', 'lfc2', 'logfold2', 'log2fc', 'logfoldchange', 'logfold', 'lf', 'logfc', 'foldchange', 'fc', 'lfc', 'enrichment']
    possible_pval_names = ['padj', 'adjp', 'pvalueadj', 'adjpvalue', 'pvaladj', 'adjpval', 'pvadj', 'adjpv', 'pvalue', 'pval', 'pv']
    
    for column, value in row.items():
        if value:  # Only process non-empty values
            column_lower = column.lower()
            if any(name in column_lower for name in possible_gene_names):
                print(f"Gene info for {value} added")
                new_triples, gene_added = process_gene(row_node, column, value, gene_added, row)
            elif any(name in column_lower for name in possible_log_names):
                new_triples, log_added = process_log_fold(row_node, column, value, log_added)
            elif any(name in column_lower for name in possible_pval_names):
                new_triples, pval_added = process_pvalue(row_node, column, value, pval_added)
            else:
                new_triples = process_other(row_node, column, value)
            triples.extend(new_triples)
    
    return triples

def process_regular_csv(csv_file_path, graph):
    dataset_uri = setup_graph(csv_file_path, graph)
    
    with open(csv_file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            triples = process_row(row, dataset_uri)
            for triple in triples:
                graph.add(triple)
    
    return graph

def process_metadata_csv(csv_file_path, graph):
    with open(csv_file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if 'pmid' in row and row['pmid']:
                pmid_uri = PMC[row['pmid']]
                for column, value in row.items():
                    if column == 'doi' and value:
                        graph.add((pmid_uri, DCT.identifier, DOI[value]))
                    elif column == 'title' and value:
                        graph.add((pmid_uri, DCT.title, Literal(value)))
                    elif column == 'year' and value:
                        graph.add((pmid_uri, DCT.date, Literal(value)))
                    elif column == 'journal' and value:
                        graph.add((pmid_uri, DCT.publisher, Literal(value)))
    return graph


if __name__ == "__main__":
    # Create the graph object which holds the triples
    graph = Graph()

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
    root_dir = "C:\\Users\\tamjh\\CodeProjects\\ASDProject\\test_data"

    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.endswith('.csv'):
                csv_file_path = os.path.join(dirpath, filename)
                print(f"Processing file: {csv_file_path}")
                if filename == 'asd_article_metadata.csv':
                    process_metadata_csv(csv_file_path, graph)
                else:
                    process_regular_csv(csv_file_path, graph)

    graph.serialize(destination='main_graph.nt', format='nt', encoding="utf-8")
    print("Combined graph has been serialized to main_graph.nt")

    #get_gene_id.cache_clear()