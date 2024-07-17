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
import urllib.request
# import textract
import pandas as pd
from functools import reduce
from metapub import PubMedFetcher
from selenium import webdriver
from urllib.parse import urljoin


# retrieve articles and convert to DOIs (which will then be converted to urls)
# Article search, returning PMIDs for articles with search terms taken from related article titles.

def get_search_result():
    Entrez.email = "thill09@student.bbk.ac.uk"
    handle = Entrez.esearch(db='pubmed',
                            term='((autism[title] or ASD[title] AND brain AND transcriptomic AND expression AND rna AND sequencing NOT Review[Publication Type]))',
                            retmax='10',
                            sort='relevance',
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
"""def get_urls(dlist):
    url_list = []
    for d in range(len(dlist)):
        prefix = 'https://doi.org/'
        new_url = prefix + dlist[d]
        url_list.append(new_url)
    return url_list"""

#https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/31097668/
def get_urls(plist):
    url_list = []
    for p in range(len(plist)):
        prefix = 'https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/'
        new_url = prefix + plist[p]
        url_list.append(new_url)
    print(url_list)
    return url_list

# retrieve supplementary files from the article
"""def get_tables(url, doi):
    main_dir = 'data'
    supp_output_dir = 'supp_data'
    new_doiref = doi.replace("/", "_")
    new_dir = new_doiref
    new_path = os.path.join(main_dir, supp_output_dir, new_dir)
    # create the directory if it doesn't exist
    os.makedirs(new_path, exist_ok=True)
    headers = {'User-Agent': 'Mozilla/5.0'}
    
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req) as u:
            html = u.read().decode('utf-8')
    except urllib.error.HTTPError as e:
        print(f"HTTP Error {e.code}: {e.reason}")
        return

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
    return"""

def get_tables(url, pmid):
    main_dir = 'data'
    supp_output_dir = 'supp_data'
    new_dir = pmid
    new_path = os.path.join(main_dir, supp_output_dir, new_dir)
    # create the directory if it doesn't exist
    os.makedirs(new_path, exist_ok=True)
    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'}
    
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req) as u:
            html = u.read().decode('utf-8')
    except urllib.error.HTTPError as e:
        print(f"HTTP Error {e.code}: {e.reason}")
        return

    soup = BeautifulSoup(html, "html.parser")
    for link in soup.find_all('a', href=True):
        href = link['href']
        if any(href.lower().endswith(x) for x in ['.csv', '.xls', '.xlsx']):
            full_url = urljoin(url, href)
            filename = os.path.join(new_path, href.rsplit('/', 1)[-1])
            if os.path.isfile(filename):
                print(f"File '{filename}' already exists. Skipping file creation.")
            else:
                print(f"Downloading {full_url} to {filename}...")
                try:
                    urllib.request.urlretrieve(full_url, filename)
                    print("Done.")
                except Exception as e:
                    print(f"Error downloading {full_url}: {e}")
    
    if not os.listdir(new_path):
        print("No files were downloaded.")
    return

def get_pdfs(url):
    main_dir = 'data'
    art_output_dir = 'article_data'
    output_path = os.path.join(main_dir, art_output_dir)
    
    # Create the directory if it doesn't exist
    os.makedirs(output_path, exist_ok=True)

    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'}
    
    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()
        html = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error fetching the page: {e}")
        return None

    soup = BeautifulSoup(html, "html.parser")
    
    # Look for PDF link in meta tags
    pdf_meta_tag = soup.find('meta', {'name': lambda name: name and 'pdf-link' in name.lower()})
    
    if pdf_meta_tag:
        pdf_url = pdf_meta_tag.get('content')
    else:
        # If not found in meta tags, look for PDF links in <a> tags
        pdf_link = soup.find('a', href=lambda href: href and href.lower().endswith('.pdf'))
        if pdf_link:
            pdf_url = pdf_link['href']
        else:
            print("No article PDF link found on the page.")
            return None

    # Ensure we have a full URL
    pdf_url = urljoin(url, pdf_url)
    
    # Extract filename
    filename_part = pdf_url.split('=')[-1] if '=' in pdf_url else pdf_url.split('/')[-1]
    if not filename_part.lower().endswith('.pdf'):
        filename_part += '.pdf'

    filename = os.path.join(output_path, filename_part)

    if os.path.isfile(filename):
        print(f"File '{filename}' already exists. Skipping download.")
        return filename

    print(f"Downloading PDF from {pdf_url}")
    try:
        response = requests.get(pdf_url, headers=headers)
        response.raise_for_status()
        with open(filename, 'wb') as fw:
            fw.write(response.content)
        print(f"PDF downloaded and saved as '{filename}'")
        return filename
    except requests.exceptions.RequestException as e:
        print(f"Error downloading PDF: {e}")
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
    url_data = get_urls(pmid_data)
    for u in url_data:
        try:
            get_pdfs(u)
        except urllib.error.HTTPError:
            pass
    for u, p in zip(url_data, pmid_data):
        try:
            get_tables(u, p)
        except urllib.error.HTTPError:
            pass
    print("All articles and data retrived")
    return 

#add a function to be used to pull only new data

if __name__ == "__main__":
    main()
