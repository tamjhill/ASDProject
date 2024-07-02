import codecs

def convert_encoding(input_file, output_file, input_encoding, output_encoding):
    """
    Converts the encoding of a file from one encoding to another.

    Args:
        input_file (str): The path to the input file.
        output_file (str): The path to the output file.
        input_encoding (str): The input encoding (e.g., 'utf-16', 'utf-8').
        output_encoding (str): The desired output encoding (e.g., 'utf-8', 'utf-16').
    """
    with codecs.open(input_file, 'r', encoding=input_encoding) as input_stream:
        content = input_stream.read()

    with codecs.open(output_file, 'w', encoding=output_encoding) as output_stream:
        output_stream.write(content)

input_file = 'result_full.nt'
output_file = 'result5.nt'
convert_encoding(input_file, output_file, 'utf-16', 'utf-8')
print(f"File encoding converted from UTF-16 to UTF-8 and saved to {output_file}")