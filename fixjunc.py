def replace_first_column(input_file, output_file):
    # Create a mapping of patterns to their replacements
    replacement_mapping = {
        "GK000031.3": "1",
        "GK000032.3": "2",
        "GK000033.3": "3",
        "CM000780.3": "4",
        "CM000781.3": "5",
        "CM000782.3": "6",
        "GK000034.3": "7",
        "CM000784.3": "8",
        "CM000785.3": "9",
        "CM000786.3": "10"
    }

    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            # Split the line by whitespace to extract the first column
            parts = line.strip().split()
            if parts:
                first_column = parts[0]
                # Check if the first column matches any of the patterns
                if first_column in replacement_mapping:
                    # If a match is found, replace it with the corresponding value
                    parts[0] = replacement_mapping[first_column]
                # Write the modified line to the output file
                f_out.write(' '.join(parts) + '\n')

if __name__ == "__main__":
    with open('junction1.txt','r') as r:
        lines = r.readlines()
        for line in lines:
            input_file_path = line.strip()  # Replace with the path of your input file
            output_file_path = line.strip()+"1"  # Replace with the path where you want the output file to be created
            replace_first_column(input_file_path, output_file_path)
            print("Replacement completed.")
