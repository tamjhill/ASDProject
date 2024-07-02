#This file searches rthe Entrez database for relevant papers, retreives their DOIs metadata, obtains a pdf 
# of the original paper, and all supporting data xlsx files

# importing libraries
from Bio import Entrez
from metapub.convert import pmid2doi
import os
from bs4 import BeautifulSoup,  SoupStrainer
import requests
# import re
from urllib.request import urlopen, urlretrieve
# import shutil
import urllib.request, urllib.error, urllib.parse
# import textract
import pandas as pd
from functools import reduce
from metapub import PubMedFetcher

# retrieve articles and convert to DOIs (which will then be converted to urls)
# Article search, returning PMIDs for articles with search terms taken from related article titles.

def get_search_result():
    Entrez.email = "thill09@student.bbk.ac.uk"
    handle = Entrez.esearch(db='pubmed',
                            term='((autism[title] or ASD[title] AND brain AND transcriptomic AND expression AND rna AND sequencing)',
                            retmax='10',
                            retmode='xml')
    search_results = Entrez.read(handle)
    return search_results


def get_pmids(search_res):
    # store as list of the PMIDs
    initial_list = search_res["IdList"]
    print(len(initial_list))
    return initial_list

# for each PMID, convert to DOI and add to new list
#NB - instead make this a dictionary to retain original PMID?
def get_dois(plist):
    doi_list = []
    for i in plist:
        try:
            doi_number = pmid2doi(i)
            if doi_number is None:
                continue
            else:
                doi_list.append(doi_number)
        except TypeError:
            continue
    print(len(doi_list))
    return doi_list

# convert each doi to a url
def get_urls(dlist):
    url_list = []
    for d in range(len(dlist)):
        prefix = 'https://doi.org/'
        new_url = prefix + dlist[d]
        url_list.append(new_url)
    return url_list


# retrieve supplementary files from the article
def get_tables(url, doi):
    main_dir = 'data'
    supp_output_dir = 'supp_data'
    new_doiref = doi.replace("/", "_")
    new_dir = new_doiref
    new_path = os.path.join(main_dir, supp_output_dir, new_dir)
    # Create the directory if it doesn't exist
    os.makedirs(new_path, exist_ok=True)
    u = urlopen(url)
    try:
        html = u.read().decode('utf-8')
    finally:
        u.close()
    soup = BeautifulSoup(html, "html.parser")
    for link in soup.select('a[href^="https://"]'):
        href = link.get('href')
        if not any(href.endswith(x) for x in ['.csv', '.xls', '.xlsx']):
            continue
        filename = os.path.join(new_path, href.rsplit('/', 1)[-1])
        if os.path.isfile(filename):
            print(f"File '{filename}' already exists. Skipping file creation.")
        else:
            print("Downloading %s to %s..." % (href, filename))
            urlretrieve(href, filename)
            print("Done.")
    return

def get_pdfs(url):
    main_dir = 'data'
    art_output_dir = 'article_data'
    u = urlopen(url)
    try:
        html = u.read().decode('utf-8')
    except urllib.error.HTTPError as e:
        if e.code in (..., 403, ...):
            pass
    finally:
        u.close()
    soup = BeautifulSoup(html, "html.parser")
    pdf_meta_tag = soup.find('meta', {'name': lambda name: name and 'pdf_url' in name.lower()})
    if pdf_meta_tag:
            # Extract the PDF URL from the content attribute of the meta tag
            pdf_url = pdf_meta_tag.get('content')
            if not pdf_url.endswith('.pdf'):
                pdf_url += '.pdf'
            print(pdf_url)
            filename = os.path.join(main_dir, art_output_dir, pdf_url.rsplit('/', 1)[-1])
            try:
                response = requests.get(pdf_url)
            except requests.exceptions.RequestException as e: 
                print(e, "error. URL not accessible.")
                pass
            if os.path.isfile(filename):
                print(f"File '{filename}' already exists. Skipping file creation.")
            else:
                with open(filename, 'wb') as fw:
                    fw.write(response.content)
            return
    return None

def get_metadata(plist, dlist):
    fetch = PubMedFetcher()
    articles = {}
    for pmid in plist:
        articles[pmid] = fetch.article_by_pmid(pmid)

    # Extract relevant information and create DataFrames
    titles = {}
    for pmid in plist:
        titles[pmid] = fetch.article_by_pmid(pmid).title
    Title = pd.DataFrame(list(titles.items()), columns=['pmid', 'title'])

#    authors = {}
#    for pmid in plist:
#        authors[pmid] = fetch.article_by_pmid(pmid).author_name
#    Author = pd.DataFrame(list(authors.items()), columns=['pmid', 'author_name'])

    dates = {}
    for pmid in plist:
        dates[pmid] = fetch.article_by_pmid(pmid).year
    Date = pd.DataFrame(list(dates.items()), columns=['pmid', 'year'])

    journals = {}
    for pmid in plist:
        journals[pmid] = fetch.article_by_pmid(pmid).journal 
    Journal = pd.DataFrame(list(journals.items()), columns=['pmid', 'journal'])

    abstracts = {}
    for pmid in plist:
        abstracts[pmid] = fetch.article_by_pmid(pmid).abstract
    Abstract = pd.DataFrame(list(abstracts.items()), columns=['pmid', 'abstract'])

    Doi = pd.DataFrame({'pmid': plist, 'doi': dlist})

    # Merge all DataFrames into a single one
    data_frames = [Title, Date, Journal, Doi, Abstract]
    df_merged = reduce(lambda  left, right: pd.merge(left, right, on=['pmid'], how='outer'), data_frames)

    # Export the merged DataFrame to a CSV file
    main_dir = 'data'
    file_path = os.path.join(main_dir, 'asd_article_metadata.csv')
    if os.path.isfile(file_path):
        print(f"File '{file_path}' already exists. Overwriting it.")
        df_merged.to_csv(file_path, index=False)
    else:
        df_merged.to_csv(file_path, index=False)
        print(f"File '{file_path}' created.")
    return None



def main():
    search_data = get_search_result()
    pmid_data = get_pmids(search_data)
    doi_data = get_dois(pmid_data)
    get_metadata(pmid_data, doi_data)
    url_data = get_urls(doi_data)
    for u in url_data:
        try:
            get_pdfs(u)
        except urllib.error.HTTPError:
            pass
    for u, d in zip(url_data, doi_data):
        try:
            get_tables(u, d)
        except urllib.error.HTTPError:
            pass
    print("All articles and data retrived")
    return 

#add a function to be used to pull only new data

if __name__ == "__main__":
    main()
