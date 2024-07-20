#currently not required - use downloaded dataset gene_ids.txt

import rdflib
from SPARQLWrapper import SPARQLWrapper, JSON
sparql = SPARQLWrapper()
    
g = rdflib.Graph()
#http://biomart.genenames.org/martservice/results

query = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>

PREFIX accesspoint: <https://biomart.genenames.org/martsemantics/hgnc_gene_config/ontology#>
PREFIX class: <biomart://biomart.genenames.org/martsemantics/hgnc_gene_config/ontology/class#>
PREFIX dataset: <biomart://biomart.genenames.org/martsemantics/hgnc_gene_config/ontology/dataset#>
PREFIX attribute: <biomart://biomart.genenames.org/martsemantics/hgnc_gene_config/ontology/attribute#>

SELECT ?hgnc_gene__hgnc_gene_id_1010 ?hgnc_gene__approved_symbol_1010 ?hgnc_gene__approved_name_1010 ?hgnc_gene__hgnc_alias_symbol__alias_symbol_108 ?hgnc_gene__hgnc_alias_name__alias_name_107 ?hgnc_gene__ensembl_gene__ensembl_gene_id_104 ?hgnc_gene__ncbi_gene__gene_id_1026 
FROM dataset:hgnc_gene_mart_2024_07_09
WHERE {
  ?x attribute:hgnc_gene__hgnc_gene_id_1010 ?hgnc_gene__hgnc_gene_id_1010 .
  ?x attribute:hgnc_gene__approved_symbol_1010 ?hgnc_gene__approved_symbol_1010 .
  ?x attribute:hgnc_gene__approved_name_1010 ?hgnc_gene__approved_name_1010 .
  ?x attribute:hgnc_gene__ensembl_gene__ensembl_gene_id_104 ?hgnc_gene__ensembl_gene__ensembl_gene_id_104 .
  ?x attribute:hgnc_gene__ncbi_gene__gene_id_1026 ?hgnc_gene__ncbi_gene__gene_id_1026
}
"""

results = g.query(query)
results.serialize(destination='all_gene_list.csv', format='csv')