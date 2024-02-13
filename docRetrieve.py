"""Retrieving the relevant documents acquired from pubmed, as well as info
directly from the artice urls
"""

def __init__(self, HF_API_TOKEN,
                   data_source_path=None,
                   data_text=None,
                   OPENAI_KEY=None) -> None:

    # Path of the file.
    self.data_source_path = data_source_path

    # File in string format directly
    self.data_text = data_text

    self.document = None
    # self.document_splited = None
    # self.embedding_model = None
    # self.embedding_type = None
    # self.OPENAI_KEY = OPENAI_KEY
    # self.HF_API_TOKEN = HF_API_TOKEN
    # self.db = None
    # self.llm = None
    # self.chain = None
    # self.repo_id = None

def get_document(self, data_source_type="TXT"):
# DS_TYPE_LIST= ["WEB", "PDF", "TXT"]
   
    data_source_type = data_source_type if data_source_type.upper() in DS_TYPE_LIST else DS_TYPE_LIST[0]
    if data_source_type == "TXT":
        if self.data_text:
            self.document = self.data_text
        elif self.data_source_path:
            loader = dl.TextLoader(self.data_source_path)
            self.document = loader.load()

    elif data_source_type == "PDF":
        if self.data_text:
            self.document = self.data_text
        elif self.data_source_path:
            loader = dl.PyPDFLoader(self.data_source_path)
            self.document = loader.load()

    elif data_source_type == "WEB":
        loader = dl.WebBaseLoader(self.data_source_path)
        self.document = loader.load()

    return self.document

# SPLIT_TYPE_LIST = ["CHARACTER", "TOKEN"]


def get_split(self, split_type="character", chunk_size=200, chunk_overlap=10):
  
  split_type = split_type.upper() if split_type.upper() 
                in SPLIT_TYPE_LIST else SPLIT_TYPE_LIST[0]
  
  if self.document:
  
      if split_type == "CHARACTER":
          text_splitter = ts.RecursiveCharacterTextSplitter(chunk_size=chunk_size, chunk_overlap=chunk_overlap)
      elif split_type == "TOKEN":
          text_splitter  = ts.TokenTextSplitter(chunk_size=chunk_size, chunk_overlap=chunk_overlap)
  
      
      # If you input a string as a document, we'll perform a split_text.
      if self.data_text:
          try:
              self.document_splited = text_splitter.split_text(text=self.document)
          except Exception as error:
              print( error)
  
      # If you upload a document, we'll do a split_documents.
      elif self.data_source_path:
          try:
              self.document_splited = text_splitter.split_documents(documents=self.document)
          except Exception as error:
              print( error)
  
  return self.document_splited


def get_embedding(self, embedding_type="HF", OPENAI_KEY=None):
  
  if not self.embedding_model:
  
    embedding_type = embedding_type.upper() if embedding_type.upper() in EMBEDDING_TYPE_LIST else EMBEDDING_TYPE_LIST[0]
  
    # If we choose to use the Hugging Face model for embedding
    if embedding_type == "HF":
        self.embedding_model = embeddings.HuggingFaceEmbeddings()

    # If we opt for the OpenAI model for embedding
    elif embedding_type == "OPENAI":
        self.OPENAI_KEY = self.OPENAI_KEY if self.OPENAI_KEY else OPENAI_KEY
        if self.OPENAI_KEY:
            self.embedding_model = embeddings.OpenAIEmbeddings(openai_api_key=OPENAI_KEY)
        else:
            print("You need to introduce a OPENAI API KEY")
  
    # The object
    self.embedding_type = embedding_type
  
    return self.embedding_model
  

  # VECTORSTORE_TYPE_LIST = ["FAISS", "CHROMA", "SVM"]

def get_storage(self, 
                 vectorstore_type = "FAISS",
                 embedding_type="HF",
                 OPENAI_KEY=None):

  self.embedding_type = self.embedding_type if self.embedding_type else embedding_type
  vectorstore_type = vectorstore_type.upper() if vectorstore_type.upper() in VECTORSTORE_TYPE_LIST else VECTORSTORE_TYPE_LIST[0]
  
  # Here we make the call to the algorithm 
  # that performed the embedding and create the object
  self.get_embedding(embedding_type=self.embedding_type, OPENAI_KEY=OPENAI_KEY)
  
  # Here we choose the type of vector store that we want to use
  if vectorstore_type == "FAISS":
      model_vectorstore = vs.FAISS
  
  elif vectorstore_type == "CHROMA":
      model_vectorstore = vs.Chroma
  
  elif vectorstore_type == "SVM":
      model_vectorstore = retrievers.SVMRetriever


  # Here we create the vector store. In this case, 
  # the document comes from raw text.
  if self.data_text:
      try:
          self.db = model_vectorstore.from_texts(self.document_splited,
                                                 self.embedding_model)
      except Exception as error:
          print( error)
  
  # Here we create the vector store. In this case,
  # the document comes from a document like pdf txt...
  elif self.data_source_path:
      try:
          self.db = model_vectorstore.from_documents(self.document_splited,
                                                     self.embedding_model)
      except Exception as error:
          print( error)
  
  return self.db


# Depending on the type of vector store we've built, 
# we'll use a specific function. All of them return
# a list of the most relevant splits.

def get_search(self, question, with_score=False):

  relevant_docs = None
  
  if self.db and "SVM" not in str(type(self.db)):
  
      if with_score:
          relevant_docs = self.db.similarity_search_with_relevance_scores(question)
      else:
          relevant_docs = self.db.similarity_search(question)

  elif self.db:
      relevant_docs = self.db.get_relevant_documents(question)
  
  return relevant_docs

    def do_question(self, 
                     question,
                     repo_id="declare-lab/flan-alpaca-large", 
                     chain_type="stuff", 
                     relevant_docs=None, 
                     with_score=False, 
                     temperature=0, 
                     max_length=300):
  # We get the most relevant splits.
  relevant_docs = self.get_search(question, with_score=with_score)
  
  # We define the LLM that we want to use, 
  # we must introduce the repo id since we are using huggingface.
  self.repo_id = self.repo_id if self.repo_id is not None else repo_id
  chain_type = chain_type.lower() if chain_type.lower() in CHAIN_TYPE_LIST else CHAIN_TYPE_LIST[0]
  
  # This check is necessary since we can call the function several times,
  # but it would not make sense to create an LLM every time the call is made.
  # So it checks if an llm already exists inside the class 
  # or if the repo_id (the type of llm) has changed.
  if (self.repo_id != repo_id ) or (self.llm is None):
     self.repo_id = repo_id 

     # We created the LLM.
     self.llm = HuggingFaceHub(repo_id=self.repo_id,huggingfacehub_api_token=self.HF_API_TOKEN,
                                model_kwargs=
                                {"temperature":temperature,
                                 "max_length": max_length})
  # We create the prompt
  prompt_template = """Use the following pieces of context to answer the question at the end. 
  If you don't know the answer, just say that you don't know, don't try to make up an answer.
  If the question is similar to [Talk me about the document], 
  the response should be a summary commenting on the most important points about the document
  
  
  {context}
  Question: {question}
  """

  PROMPT = PromptTemplate(
    template=prompt_template, input_variables=["context", "question"]
  )  
  
  # We create the chain, chain_type= "stuff".
  self.chain = self.chain if self.chain is not None 
                else load_qa_chain(self.llm,
                                   chain_type=chain_type,
                                   prompt = PROMPT)
  
  # We make the query to the LLM using the prompt
  # We check if there is a chain already defined, 
  # if it does not exist it is created
  response = self.chain({"input_documents": relevant_docs, "question": question}, return_only_outputs=True)
  
  return response