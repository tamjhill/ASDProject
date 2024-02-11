
# importing libraries
from Bio import Entrez
from metapub.convert import pmid2doi
import os
from bs4 import BeautifulSoup,  SoupStrainer
import requests
# import re
from urllib.request import urlopen, urlretrieve
import urllib.parse
# import shutil
import os
import urllib

# retrieve articles and convert to DOIs (which will then be converted to urls)
# Article search, returning PMIDs for articles with search terms taken from related article titles.


def get_articles():
    Entrez.email = "thill09@student.bbk.ac.uk"
    handle = Entrez.esearch(db='pubmed',
                            term='(autism[title] AND brain AND transcriptomic AND expression AND rna)',
                            retmax='500',
                            retmode='xml')
    search_results = Entrez.read(handle)

    # store as list of the PMIDs
    initial_list = search_results["IdList"]
    print(len(initial_list))
    return initial_list

# for each PMID, convert to DOI and add to new list

#NB - instead make this a dictionary to retain original PMID
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
    testurl = url_list[2]
    return testurl


def find_data(url):
    output_dir = 'data'
    u = urlopen(url)
    try:
        html = u.read().decode('utf-8')
    finally:
        u.close()
    soup = BeautifulSoup(html, "html.parser")
    for link in soup.select('a[href^="https://"]'):
        href = link.get('href')
        if not any(href.endswith(x) for x in ['.csv', '.xls', '.xlsx', '.txt']):
            continue
        filename = os.path.join(output_dir, href.rsplit('/', 1)[-1])
        print("Downloading %s to %s..." % (href, filename))
        urlretrieve(href, filename)
        print("Done.")
    return

find_data('https://doi.org/10.1038/s41586-023-06473-y')


def main():
    art_data = get_articles()
    doi_data = get_dois(art_data)
    url_data = get_urls(doi_data)
    find_data(url_data)


if __name__ == "__main__":
    main()

#https://doi.org/10.1038/s41586-023-06473-y

# t_url = 'https://doi.org/10.1016/j.molcel.2016.11.033'
#t2_url = 'https://doi.org/10.1038/s41586-023-06473-y'
#find_data(t2_url)