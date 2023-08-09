def read_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    return [line.strip() for line in lines]  # Strip leading/trailing whitespace from each line

def write_file(filename, lines):
    with open(filename, 'w') as file:
        file.writelines(lines)

def grep_and_rename_columns(data_lines, rename_lines):
    # Create a dictionary to store the rename mappings from the rename file
    rename_dict = {}
    for line in rename_lines:
        old_column, new_column = line.strip().split()
        rename_dict[old_column] = new_column

    # Extract the specific columns from the data
    header_line = data_lines[0].split()
    data_lines = data_lines[1:]  # Remove the first line (header) from the data
    grep_column_indexes = [i for i in range(len(header_line)) if header_line[i] in rename_dict]

    # Process the data and update the column names
    updated_header_line = []
    for col in header_line:
        if col in rename_dict:
            updated_header_line.append(rename_dict[col])
    output_lines = [' '.join(updated_header_line) + '\n']  # Include the updated header line in the output

    for line in data_lines:
        values = line.split()
        output_line_values = [values[i] for i in grep_column_indexes]
        output_line = ' '.join(output_line_values) + '\n'
        output_lines.append(output_line)

    return output_lines

# Example usage:
data_file = 'all001_perind.counts'
rename_file = 'rtSRR_lines.txt'
output_file = 'output.txt'

# Read data from files A and B
data_lines = read_file(data_file)
rename_lines = read_file(rename_file)

# Get the updated header line with the replaced column names and the corresponding data values
updated_lines = grep_and_rename_columns(data_lines, rename_lines)

# Write the updated header line and data to the output file
write_file(output_file, updated_lines)
print("Header line with updated column names and data written to output file.")
