import argparse
import ast
from collections import OrderedDict
import itertools
import sys

parser = argparse.ArgumentParser(description="""
Python script parses through BLAST's tabular output, and filters out hits that
are below a cutoff e-value (thus removing the need to re-run searches with
varying e-values).

The output of the script is XML (--xml), table & each hit on one line (--table),
table & each hsp on each table (--vtable).""")

parser.add_argument('blast_tabular', metavar="tabular_filename",
                    type=argparse.FileType('r'), nargs='?',
                    default=sys.stdin, help="BLAST tabular filename.")
parser.add_argument('-e', '--e_value', metavar='e_value', type=float,
                    help='set maximum E value.')
parser.add_argument('-b', '--bit_score', metavar='bit_score', type=float,
                    help='set minimum total bit score value.')
parser.add_argument('-c', '--coverage', metavar='coverage_pct', type=float,
                    help='set minimum coverage (in %%) for individual hits.')
parser.add_argument('-i', '--identity', metavar='identity_pct', type=float,
                    help='set minimum local identity (in %%).')
parser.add_argument('-t', '--top', metavar='n', type=int,
                    help='select for top n hits.')
parser.add_argument('-C', '--compress', action='store_true',
                    help="""print compressed output, showing all hits for each
                    query on a single line. Can be _very_ verbose.""")
parser.add_argument('-N', '--remove_N', metavar='query_file',
                    type=argparse.FileType('r'),
                    help="""discards Ns before calculating coverage %% (best 
                    for scaffolds). Requires FASTA file of query sequences.""")
fasta_opt = parser.add_mutually_exclusive_group(required=True)
fasta_opt.add_argument('--xml', action='store_const', dest='output_format',
                       const='xml', help='produce XML output.')
fasta_opt.add_argument('--table', action='store_const', dest='output_format',
                       const='table', help='produce tabular output.')
fasta_opt.add_argument('--vtable', action='store_const', dest='output_format',
                       const='vtable',
                       help='produce verbose tabular output (1 line per Hsp).')
parser.add_argument('--noheader', action='store_true',
                    help='disable printing of header.')
args = parser.parse_args()

def compress_output(output):
    c_output = ''
    if not args.noheader:
        c_output += '\t'.join(['Query', 'Hit accession', 'Hit description',
                               'Query length', 'Hit length', 
                               'Query (start, end)', 'Hit (start, end)',
                               'Frame', 'Max bit score', 'Total bit score',
                               'Identity', 'Identity %',
                               'Coverage %', 'Expect']) + '\n'
    
    # create a dictionary to store lines in the output that belong to each query
    lines_per_query = OrderedDict()
    
    for line in output.split('\n')[1:]:
        if not line:
            break
        
        cols = line.split('\t')
        if cols[0] not in lines_per_query:
            lines_per_query[cols[0]] = []
        
        lines_per_query[cols[0]].append(line)
    
    return lines_per_query

# Load the tabular data
with open(args.blast_tabular.name, 'r') as file:
    tabular_data = file.read()

# Process the tabular data
processed_data = compress_output(tabular_data)

# Output the processed data based on the selected output format
if args.output_format == 'table':
    for query, lines in processed_data.items():
        for line in lines:
            print(line)
elif args.output_format == 'vtable':
    for query, lines in processed_data.items():
        for line in lines:
            print(line)
elif args.output_format == 'xml':
    # Convert processed_data to XML format and output
    pass  # Placeholder, actual implementation required
