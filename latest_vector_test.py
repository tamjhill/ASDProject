from langchain_community.document_loaders import TextLoader
from langchain_community.document_loaders import UnstructuredFileLoader
from langchain.text_splitter import CharacterTextSplitter
from langchain_community.vectorstores import FAISS
from langchain_community.vectorstores import FAISS
from langchain_community.embeddings import OpenAIEmbeddings
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.document_loaders.sitemap import SitemapLoader
from langchain.embeddings import HuggingFaceEmbeddings
from sentence_transformers import SentenceTransformer, util
from tika import parser




filepath = 'C:/Users/tamjh/CodeProjects/ASDProject/article_data/s12859-023-05278-0.pdf'
parsed_document = parser.from_file(filepath)
print(parsed_document['content'])
#You can also see all available attributes by using the following line:

print(parsed_document.keys())


#better loader for pdfs? nougat?
#loader = UnstructuredFileLoader(
#    "layout-parser-paper-fast.pdf", strategy="hi_res", mode="elements"
#)
#
#docs = loader.load()