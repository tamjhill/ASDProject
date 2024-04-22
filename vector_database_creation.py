import os
from pypdf import PdfReader
import openpyxl
from sentence_transformers import SentenceTransformer
import faiss

# Initialize the sentence transformer model (to generate embeddings)
model = SentenceTransformer('all-MiniLM-L6-v2')

# Initialize the FAISS index (to create the vector database)
index = faiss.IndexFlatIP(model.get_sentence_embedding_dimension())
index_ids = []
index_metadata = []

# Iterate over files in a directory (inc. sub-directories)
data_dir = 'data'
for root, dirs, files in os.walk(data_dir):
    for filename in files:
        file_path = os.path.join(root, filename)
        
        if filename.endswith('.pdf'):
            # Process PDF files
            reader = PdfReader(file_path)
            text = ''
            for page in reader.pages:
                text += page.extract_text()
            embedding = model.encode([text])
            index.add(embedding)
            index_ids.append(len(index_ids))
            index_metadata.append({'file': filename, 'type': 'pdf', 'directory': root})
            
        elif filename.endswith('.xlsx'):
            # Process Excel files
            workbook = openpyxl.load_workbook(file_path)
            text = ''
            for sheet in workbook.worksheets:
                for row in sheet.iter_rows():
                    for cell in row:
                        text += str(cell.value) + ' '
            embedding = model.encode([text])
            index.add(embedding)
            index_ids.append(len(index_ids))
            index_metadata.append({'file': filename, 'type': 'xlsx', 'directory': root})
            
        elif filename.endswith('.txt'):
            # Process text files
            with open(file_path, 'r') as f:
                text = f.read()
            embedding = model.encode([text])
            index.add(embedding)
            index_ids.append(len(index_ids))
            index_metadata.append({'file': filename, 'type': 'txt', 'directory': root})

# Save the FAISS index
faiss.write_index(index, 'index.faiss')

# Save the index metadata
import pickle
with open('metadata.pkl', 'wb') as f:
    pickle.dump(index_metadata, f)