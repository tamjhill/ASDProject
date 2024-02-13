# import os
# import torch
from transformers import (
    AutoTokenizer,
    AutoModelForCausalLM,
    pipeline,
    AutoConfig
)

from langchain.text_splitter import CharacterTextSplitter
from langchain_community.document_transformers import Html2TextTransformer
from langchain_community.document_loaders import AsyncChromiumLoader
from langchain.embeddings.huggingface import HuggingFaceEmbeddings
from langchain_community.vectorstores import FAISS

from langchain.prompts import PromptTemplate
from langchain.schema.runnable import RunnablePassthrough
from langchain_community.llms import HuggingFacePipeline
from langchain.chains import LLMChain

import nest_asyncio

#################################################################
# Tokenizer
#################################################################

model_name = 'mistralai/Mistral-7B-Instruct-v0.1'

model_config = AutoConfig.from_pretrained(
    model_name,
)

tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
tokenizer.pad_token = tokenizer.eos_token
tokenizer.padding_side = "right"

#################################################################
# Load pre-trained config
#################################################################
model = AutoModelForCausalLM.from_pretrained(
    model_name
)

text_generation_pipeline = pipeline(
    model=model,
    tokenizer=tokenizer,
    task="text-generation",
    temperature=0,
    repetition_penalty=1.1,
    return_full_text=True,
    max_new_tokens=1000,
)

mistral_llm = HuggingFacePipeline(pipeline=text_generation_pipeline)

# Creating the vector database for the RAG to search:

nest_asyncio.apply()

# Articles to index
articles = ["https://doi.org/10.1038/s41586-023-06473-y"]

# Scrapes the article above
loader = AsyncChromiumLoader(articles)
docs = loader.load()

# Converts HTML to plain text
html2text = Html2TextTransformer()
docs_transformed = html2text.transform_documents(docs)

# Chunk text
text_splitter = CharacterTextSplitter(chunk_size=500,
                                      chunk_overlap=50)
chunked_documents = text_splitter.split_documents(docs_transformed)

# Load chunked documents into the FAISS index
db = FAISS.from_documents(chunked_documents,
                          HuggingFaceEmbeddings(model_name='sentence-transformers/all-mpnet-base-v2'))

# Link between the vector db, Langchain and the LLM
retriever = db.as_retriever(
    search_type="similarity",
    search_kwargs={'k': 6}
)

# Create prompt template - context provided by the chunks taken from the vector DB
prompt_template = """
### [INST] Instruction: Answer the questions on this article about Autism Spectrum Disorder using the following context:

{context}

### QUESTION:
{question} [/INST]
 """

# Create prompt from prompt template
prompt = PromptTemplate(
    input_variables=["context", "question"],
    template=prompt_template,
)

# Create llm chain
llm_chain = LLMChain(llm=mistral_llm, prompt=prompt)

rag_chain = (
        {"context": retriever, "question": RunnablePassthrough()}
        | llm_chain
)

rag_chain.invoke("Give a summary of this article and state the key genes/transcripts discussed and their roles")





