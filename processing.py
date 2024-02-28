
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


#from_url from module 'pdfkit' 

"""
url = "https://doi.org/10.1038/s41586-023-06473-y"
my_dir = "data" 
my_file = "testfile1.pdf"
new_file = os.path.join(my_dir, my_file)
response = requests.get(url)
with open('testfile1.pdf', 'wb') as f:
    f.write(response.content)
    print("file created")"""

def get_search_result():
    Entrez.email = "thill09@student.bbk.ac.uk"
    handle = Entrez.esearch(db='pubmed',
                            term='(autism[title] AND brain AND transcriptomic AND expression AND rna)',
                            retmax='500',
                            retmode='xml')
    search_results = Entrez.read(handle)
    return search_results

#def get_pdfs(search_res):
#    pdf_list = []
#    urlretrieve(url, any_path)

#with open(another_path, "w") as textfile:
#    textfile.write(textract.process(
#        any_path,
#        extension='pdf',
#        method='pdftotext',
#        encoding="utf_8",
#    ))
    #for record in search_res:
     #   if record.get('MedlineCitation'):
      ##      if record['MedlineCitation'].get('OtherID'):
        #        for other_id in record['MedlineCitation']['OtherID']:
         #           if other_id.title().startswith('Pmc'):
          #              new_title = other_id.title().upper()
           #             pdf_link = (f'http://www.ncbi.nlm.nih.gov/pmc/articles/{new_title}/pdf/')
            #            pdf_list.append(pdf_link)
             #       else:
              #          continue
#    print(len(pdf_list))
 #   return pdf_list

def get_pmids(search_res):
    # store as list of the PMIDs
    initial_list = search_res["IdList"]
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


# retrieve supplementary files from the article
def get_tables(url):
    supp_output_dir = 'supp_data'
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
        filename = os.path.join(supp_output_dir, href.rsplit('/', 1)[-1])
        print("Downloading %s to %s..." % (href, filename))
        urlretrieve(href, filename)
        print("Done.")
    return

def get_pdfs(url):
    art_output_dir = 'article_data'
    u = urlopen(url)
    try:
        html = u.read().decode('utf-8')
    finally:
        u.close()
    soup = BeautifulSoup(html, "html.parser")
    for pdf in soup.select('meta pdf'):
        if pdf.get('meta pdf') == None:
            print("none")
            continue
        my_meta = pdf.get('my_meta')
        if not any(my_meta.endswith(x) for x in ['.pdf']):
            print("no pdf")
            continue
        filename = os.path.join(art_output_dir, (my_meta.rsplit('/', 1)[-1]))
        #print("Downloading %s to %s..." % (my_meta, filename))
        urlretrieve(my_meta, filename)
        print("Done.")
    return

def get_pdf_url(url):
    print("hello")
    # Send a GET request to the URL
    response = requests.get(url)
    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        # Parse the HTML content with BeautifulSoup
        soup = BeautifulSoup(response.text, 'html.parser')

        # Find the relevant meta tag with the PDF URL
        pdf_meta_tag = soup.find('meta', {'name': lambda name: name and 'pdf_url' in name.lower()})

        if pdf_meta_tag:
            # Extract the PDF URL from the content attribute of the meta tag
            pdf_url = pdf_meta_tag.get('content')

            # You may want to validate or manipulate the URL if necessary
            # ...
            print(pdf_url)
            return pdf_url

    # If the request was not successful, print an error message
    print(f"Error: Unable to fetch content from {url}")
    return None


"""
        # Find the relevant meta tag or other HTML element containing the PDF URL
        # Adjust the following line based on the structure of the HTML
        pdf_url = soup.find('meta', {'name': '*pdf_url'})
        if pdf_url:
            # Extract the PDF URL from the attribute (e.g., content attribute for meta tag)
            pdf_url = pdf_url.get('content')
            filename = os.path.join(art_output_dir, href.rsplit('/', 1)[-1])
            #print("Downloading %s to %s..." % (href, filename))
            urlretrieve(pdf_url, filename)
            # You may want to validate or manipulate the URL if necessary
            # ...

    return
       """ 

"""
title = soup.find("meta", property="citation_pdf_url")
#url = soup.find("meta", property="og:url")

print(title["content"] if title else "No meta title given")
print(url["content"] if url else "No meta url given")
        # Get response object for link
        #response = requests.get(link.get('href'))
    return

# find_data('https://doi.org/10.1038/s41586-023-06473-y')"""


def main():
    #search_data = get_search_result()
    #supp_data = get_tables('https://doi.org/10.1038/s41586-023-06473-y')
    #pdf_data = get_pdfs('https://doi.org/10.1038/s41586-023-06473-y')
    #pmid_data = get_pmids(search_data)
    #doi_data = get_dois(art_data)
    #url_data = get_urls(doi_data)
    #find_data(url_data)
    get_pdf_url('https://doi.org/10.1038/s41586-023-06473-y')
    return 


if __name__ == "__main__":
    main()

#https://doi.org/10.1038/s41586-023-06473-y

# t_url = 'https://doi.org/10.1016/j.molcel.2016.11.033'
#t2_url = 'https://doi.org/10.1038/s41586-023-06473-y'
#find_data(t2_url)-