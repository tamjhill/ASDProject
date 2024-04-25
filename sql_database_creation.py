
"""creating seperate sql db fpor specific log fold data for extra analysis"""
#checking which files have been stored:
# Load the metadata
with open('metadata.pkl', 'rb') as f:
    index_metadata = pickle.load(f)

# Print the metadata for each file
#for metadata in index_metadata:
#    print(f"File: {metadata['directory']}/{metadata['file']} (Type: {metadata['type']})")

# Load the FAISS index and metadata
index = faiss.read_index('index.faiss')
with open('metadata.pkl', 'rb') as f:
    index_metadata = pickle.load(f)

# Initialize the sentence transformer model
model = SentenceTransformer('all-MiniLM-L6-v2')

# Perform a similarity search for "finance report"
query = "fold change"
query_vector = model.encode([query])
k = 10  # Number of nearest neighbors to retrieve
distances, indices = index.search(query_vector, k)

# Extract the desired information from the relevant files
fold_change_data = []
for index_id in indices[0]:
    metadata = index_metadata[index_id]
    file_path = os.path.join(metadata['directory'], metadata['file'])
    
    # Load and parse the file to extract the desired information
    report_data = extract_finance_report(file_path)
    fold_change_data.append(report_data)

# Create a new SQLite database and store the extracted information
conn = sqlite3.connect('finance_reports.db')
c = conn.cursor()
c.execute('''CREATE TABLE reports
             (report_id INTEGER PRIMARY KEY, report_date TEXT, revenue REAL, expenses REAL)''')

for report in finance_reports:
    c.execute("INSERT INTO reports (report_date, revenue, expenses) VALUES (?, ?, ?)",
              (report['date'], report['revenue'], report['expenses']))

conn.commit()
conn.close()