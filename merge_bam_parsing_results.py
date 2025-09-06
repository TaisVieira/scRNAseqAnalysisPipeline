#!/usr/bin/env python3
import argparse
import os
from collections import defaultdict

def merge_unique_reads_files(input_files, output_file):
    """
    Merge multiple unique_reads_per_cell.txt files into a single output file.
    
    Args:
        input_files (list): List of input file paths
        output_file (str): Output file path
    """
    print(f"Merging {len(input_files)} unique reads files into {output_file}")
    
    # Dictionary to store merged data: barcode -> total read count
    merged_data = defaultdict(int)
    
    # Total number of records processed
    total_records = 0
    
    # Process each input file
    for file_path in input_files:
        print(f"Processing {file_path}")
        
        with open(file_path, 'r') as f:
            # Skip header line
            header = f.readline().strip()
            
            # Process each line
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) != 2:
                    continue
                    
                barcode, read_count = parts
                try:
                    merged_data[barcode] += int(read_count)
                    total_records += 1
                except ValueError:
                    print(f"Warning: Could not parse read count from line: {line}")
    
    # Write merged data to output file
    with open(output_file, 'w') as out:
        print("Barcode\tReads", file=out)
        
        for barcode, count in merged_data.items():
            print(f"{barcode}\t{count}", file=out)
    
    print(f"Merged {total_records} records from {len(input_files)} files")
    print(f"Output contains {len(merged_data)} unique barcodes")

def merge_duplicate_reads_files(input_files, output_file):
    """
    Merge multiple duplicate_reads.txt files into a single output file.
    
    Args:
        input_files (list): List of input file paths
        output_file (str): Output file path
    """
    print(f"Merging {len(input_files)} duplicate reads files into {output_file}")
    
    # Dictionary to store merged data: sequence -> duplicate count
    merged_data = defaultdict(int)
    
    # Total number of records processed
    total_records = 0
    
    # Process each input file
    for file_path in input_files:
        print(f"Processing {file_path}")
        
        with open(file_path, 'r') as f:
            # Skip header line
            header = f.readline().strip()
            
            # Process each line
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) != 2:
                    continue
                    
                sequence, duplicate_count = parts
                try:
                    merged_data[sequence] += int(duplicate_count)
                    total_records += 1
                except ValueError:
                    print(f"Warning: Could not parse duplicate count from line: {line}")
    
    # Write merged data to output file
    with open(output_file, 'w') as out:
        print("Mapped_sequences\tDuplicates", file=out)
        
        for sequence, count in merged_data.items():
            print(f"{sequence}\t{count}", file=out)
    
    print(f"Merged {total_records} records from {len(input_files)} files")
    print(f"Output contains {len(merged_data)} unique mapped sequences")

def main():
    parser = argparse.ArgumentParser(description='Merge multiple output files from BAM processing')
    parser.add_argument('--unique-reads', '-u', nargs='+', required=True, 
                        help='List of unique_reads_per_cell.txt files to merge')
    parser.add_argument('--duplicate-reads', '-d', nargs='+', required=True,
                        help='List of duplicate_reads.txt files to merge')
    parser.add_argument('--output-prefix', '-o', default='merged',
                        help='Prefix for output files (default: merged)')
    
    args = parser.parse_args()
    
    # Create output filenames
    unique_output = f"{args.output_prefix}_unique_reads_per_cell.txt"
    duplicate_output = f"{args.output_prefix}_duplicate_reads.txt"
    
    # Verify input files exist
    for file_list in [args.unique_reads, args.duplicate_reads]:
        for file_path in file_list:
            if not os.path.exists(file_path):
                print(f"Error: Input file does not exist: {file_path}")
                return 1
    
    # Merge unique reads files
    merge_unique_reads_files(args.unique_reads, unique_output)
    
    # Merge duplicate reads files
    merge_duplicate_reads_files(args.duplicate_reads, duplicate_output)
    
    print(f"Merged results written to:")
    print(f"  {unique_output}")
    print(f"  {duplicate_output}")
    
    return 0

if __name__ == "__main__":
    exit(main())
