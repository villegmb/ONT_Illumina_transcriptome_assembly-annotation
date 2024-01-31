#!/usr/bin/env python3

import sys  # Import the sys module for writing to stderr

def read_transcripts(filename):
    transcripts = set()
    with open(filename, 'r') as file:
        for line in file:
            transcript = line.strip()
            transcripts.add(transcript)
    return transcripts

def process_fasta(input_fasta, output_fasta, transcripts_file):
    transcripts = read_transcripts(transcripts_file)
    transcript_sequences = {}

    with open(input_fasta, 'r') as input_file, open(output_fasta, 'w') as output_file:
        name = None
        sequence = []

        for line in input_file:
            line = line.strip()

            if line.startswith('>'):
                # Store previous sequence (if any)
                if name is not None and ''.join(sequence) != '':
                    transcript = name.split('.')[0:3]
                    print("transcript",transcript, file=sys.stderr)
                    transcript = '.'.join(transcript)
                    if transcript in transcripts:
                        output_file.write(f'>{name}\n')
                        output_file.write(''.join(sequence) + '\n')
                        print("sequence",sequence, file=sys.stderr)
                # Reset variables for new sequence
                name = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)

        # Write the last sequence
        if name is not None and ''.join(sequence) != '':
            transcript = name.split('.')[0:3]
            transcript = '.'.join(transcript)
            if transcript in transcripts:
                output_file.write(f'>{name}\n')
                output_file.write(''.join(sequence) + '\n')

# File paths

input_fasta = "/ibex/tmp/c2078/Heat_stress_analysis/scripts/Ahem_transcripts_with_reference.fa.transdecoder_dir/longest_orfs.pep"
transcripts_file = "/ibex/tmp/c2078/Heat_stress_analysis/scripts/output_column_2_no_quotes.txt"
output_fasta = "/ibex/tmp/c2078/Heat_stress_analysis/scripts/Ahem_transcripts_with_reference.fa.transdecoder_dir/longest_orfs.pep_filtered_protein.f
asta"

# Process the fasta file
process_fasta(input_fasta, output_fasta, transcripts_file)
