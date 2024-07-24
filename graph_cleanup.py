import re

def process_nt_file(input_file, output_file):
    # regex's to find decimal numbers as strings, and empty literals
    decimal_pattern = re.compile(r'"(-?\d+(\.\d+)?([eE][-+]?\d+)?)"')
    empty_literal_pattern = re.compile(r'\s+"(\s*|-)"(\^\^<[^>]+>)?\s*\.$')

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Skip empty lines or lines with empty literals
            if line.strip() == '' or empty_literal_pattern.search(line):
                continue
            
            # Check if the line contains decimal number as a string, replace with XSD double label
            match = decimal_pattern.search(line)
            if match:
                modified_line = decimal_pattern.sub(r'"\1"^^<http://www.w3.org/2001/XMLSchema#double>', line)
                outfile.write(modified_line)
            else:
                outfile.write(line)

input_file = 'main_graph.nt'
output_file = 'cleaned_maingraph.nt'
process_nt_file(input_file, output_file)