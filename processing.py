
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
import textract

# retrieve articles and convert to DOIs (which will then be converted to urls)
# Article search, returning PMIDs for articles with search terms taken from related article titles.

def get_search_result():
    Entrez.email = "thill09@student.bbk.ac.uk"
    handle = Entrez.esearch(db='pubmed',
                            term='(autism[title] AND brain AND transcriptomic AND expression AND rna)',
                            retmax='5',
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
        #new_url = prefix.replace("%3A", ":")
        url_list.append(new_url)
    # testurl = url_list[2]
    return url_list


# retrieve supplementary files from the article
def get_tables(url):
    supp_output_dir = 'supp_data'
    new_dir = url.rsplit('/', 1)[-1]
    new_path = os.path.join(supp_output_dir, new_dir)
    os.mkdir(new_path)
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
        print("Downloading %s to %s..." % (href, filename))
        urlretrieve(href, filename)
        print("Done.")
    return

def get_pdfs(url):
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
            filename = os.path.join(art_output_dir, pdf_url.rsplit('/', 1)[-1])
            try:
                response = requests.get(pdf_url)
            except requests.exceptions.RequestException as e: 
                print(e, "error. URL not accessible.")
                pass
            with open(filename, 'wb') as fw:
                fw.write(response.content)
            return
    return None


def main():
    search_data = get_search_result()
    pmid_data = get_pmids(search_data)
    doi_data = get_dois(pmid_data)
    url_data = get_urls(doi_data)
    for u in url_data:
        try:
            get_pdfs(u)
            get_tables(u)
        except urllib.error.HTTPError:
            pass
    print("All articles and data retrived")
    return 


#add a function to be used to pull only new data

if __name__ == "__main__":
    main()

#test article - https://doi.org/10.1038/s41586-023-06473-y
