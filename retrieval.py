import os
import fitz  # PyMuPDF library
import openpyxl
from sentence_transformers import SentenceTransformer
import faiss
import pickle
import re

"""returning info"""

# Load the FAISS index and metadata
index = faiss.read_index('index.faiss')
with open('metadata.pkl', 'rb') as f:
    index_metadata = pickle.load(f)

# Initialize the sentence transformer model
model = SentenceTransformer('all-MiniLM-L6-v2')

# Function to retrieve data based on a prompt
def retrieve_data(prompt, num_results=5):
    # Encode the prompt into a vector
    query_vector = model.encode([prompt])

    # Perform a similarity search
    distances, indices = index.search(query_vector, num_results)

    # Retrieve the relevant data
    relevant_data = []
    for distance, index_id in zip(distances[0], indices[0]):
        metadata = index_metadata[index_id]
        file_path = os.path.join(metadata['directory'], metadata['file'])

        # Load and process the file based on its type
        if metadata['type'] == 'pdf':
            text = process_pdf(file_path)
        elif metadata['type'] == 'xlsx':
            text = process_excel(file_path)
        elif metadata['type'] == 'txt':
            with open(file_path, 'r') as f:
                text = f.read()

        relevant_data.append({
            'text': text,
            'file': metadata['file'],
            'type': metadata['type'],
            'distance': distance
        })

    return relevant_data

# Helper functions to process different file types
def process_pdf(file_path):
    import fitz
    pdf_doc = fitz.open(file_path)
    text = ''
    for page in pdf_doc:
        text += page.get_text()
    return text

def process_excel(file_path):
    import openpyxl
    workbook = openpyxl.load_workbook(file_path)
    text = ''
    for sheet in workbook.worksheets:
        for row in sheet.iter_rows():
            for cell in row:
                text += str(cell.value) + ' '
    return text

# Example usage
prompt = "which genes are upregulated in ASD?"
relevant_data = retrieve_data(prompt)

for data in relevant_data:
    print(f"File: {data['file']} (Type: {data['type']}, Distance: {data['distance']:.4f})")
    print(data['text'][:100], "...")
    print("-" * 50)