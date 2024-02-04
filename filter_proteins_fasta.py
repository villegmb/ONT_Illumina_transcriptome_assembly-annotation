#!/usr/bin/env python3

import sys

def read_transcripts(transcripts_file):
    transcripts = set()
    with open(transcripts_file, 'r') as file:
        for line in file:
            transcripts.add(line.strip())
            print("line:", line)
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
                # Extract the identifier from the header line
                #identifier = line.split()[0][1:].split('.')[0:3]
                #identifier = '.'.join(identifier)
                identifier = line.split()[0][1:]
                print("Identifier:", identifier)

                # Check if the identifier matches any transcript name
                if identifier in transcripts:
                    if name is not None and ''.join(sequence) != '':
                        print("Name:", name)
                        print("Sequence:", ''.join(sequence))
                        print("Transcripts:", transcripts)
                        output_file.write(f'{name}\n')
                        output_file.write(''.join(sequence) + '\n')

                    name = line
                    sequence = []
                else:
                    name = None
            else:
                if name is not None:
                    sequence.append(line)

        # Write the last sequence
        if name is not None and ''.join(sequence) != '':
            print("Saving in new fasta:")
            print("Name:", name)
            print("Sequence:", ''.join(sequence))
            print("Transcripts:", transcripts)
            output_file.write(f'{name}\n')
            output_file.write(''.join(sequence) + '\n')

# Define file paths
input_fasta="/ibex/tmp/c2078/Heat_stress_analysis/scripts/Ahem_transcripts_with_reference.fa.transdecoder_dir/longest_orfs.pep"
transcripts_file="/ibex/tmp/c2078/Heat_stress_analysis/scripts/extracted_protein_names_without_GOTERM.txt"
output_fasta="/ibex/tmp/c2078/Heat_stress_analysis/scripts/Ahem_transcripts_with_reference.fa.transdecoder_dir/filtered_noGOTERM_long_proteins.fasta"

process_fasta(input_fasta, output_fasta, transcripts_file)

