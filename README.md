Outline of the project and code:

1. finding relevant articles
    - processing.py

2. retrieving article metadata, abstracts and supporting files
    - processing.py

3. checking each data file for relevant expression info
    - dataconvert.py
   
4. data cleaning
    - csv_data_cleaning.py 

5. mapping to rdf triples
    - create_rdf_triples.py
    - graph_cleanup.py

6. data testing and analysis (see analysis directory)
    - general_tests.ipynb
    - graphanalysis.ipynb
    - usecase1.ipynb
    - usecase2.ipynb
    - usecase2.ipynb

Retrieved article outputs are stored in the 'data' directory, the rdf graph is stored in 'cleaned_maingraph.nt' .
