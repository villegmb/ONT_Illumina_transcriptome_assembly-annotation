dividing_files.py
import os

# Function to split the TSV file
def split_tsv(input_file, output_folder, num_files):
    with open(input_file, 'r') as f:
        header = next(f)  # Read the header
        file_handles = []

        # Create output files and store their handles
        for i in range(num_files):
            filename = os.path.join(output_folder, f'output_{i}.tsv')
            file_handles.append(open(filename, 'w'))
            file_handles[i].write(header)  # Write header to each file

        # Distribute rows across files
        for i, line in enumerate(f):
            file_handles[i % num_files].write(line)

        # Close all files
        for handle in file_handles:
            handle.close()

if __name__ == "__main__":
    input_file = "/ibex/tmp/c2078/Heat_stress_analysis/scripts/go_analysis/ahem_vs_sprot.t20.tsv"  # Specify your input file
    output_folder = "/ibex/tmp/c2078/Heat_stress_analysis/scripts/go_analysis/ahem_vs_sprot.t20_divided"  # Specify the output folder
    num_files = 100  # Number of files to split into

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    split_tsv(input_file, output_folder, num_files)
