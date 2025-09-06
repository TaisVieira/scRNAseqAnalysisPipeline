#!/usr/bin/env python3
import argparse
import sys
import re

def extract_gff_info(input_file, output_file):
    """
    Extract specific information from a GFF file and output in simplified format.
    
    Output format:
    sequence_id gene_id start end strand
    """
    with open(input_file, 'r') as gff_file, open(output_file, 'w+') as annotation_file:
        # Skip header lines (those starting with #)
        for line in gff_file:
            if line.startswith('#'):
                continue
                
            try:
                # Process data lines
                fields = line.strip().split('\t')
                
                # Check if this is a protein_coding_gene line
                if len(fields) >= 9 and fields[2] == "protein_coding_gene":
                    seq_id = fields[0]  # Chromosome/sequence ID
                    
                    # Extract gene ID from the attributes field
                    attributes = fields[8].split(';')
                    gene_id = ""
                    for attr in attributes:
                        if attr.startswith("ID="):
                            gene_id = attr[3:]  # Remove "ID=" prefix
                            break
                    
                    # Process strand and positions with the +500bp adjustment
                    if fields[6] == "+":  # Forward strand
                        position = "foward_strand"
                        start = int(fields[3])
                        end = int(fields[4]) + 500  # Add 500bp downstream
                    elif fields[6] == "-":  # Reverse strand
                        position = "reverse_strand"
                        # For reverse strand, we add 500bp to start (which is actually 3' end genomically)
                        start = int(fields[3]) - 500  # Add 500bp downstream (5' direction in genome)
                        if start < 1:  # Ensure we don't go below position 1
                            start = 1
                        end = int(fields[4])
                    else:
                        continue  # Skip if strand is neither + nor -
                    
                    # Always ensure start < end in the output
                    if start > end:
                        start, end = end, start
                        
                    # Write to output file
                    print(f"{seq_id}\t{gene_id}\t{start}\t{end}\t{position}", file=annotation_file)
                    
            except Exception as e:
                # Skip problematic lines
                continue

def main():
    parser = argparse.ArgumentParser(description='Extract specific information from GFF file')
    parser.add_argument('input', help='Input GFF file')
    parser.add_argument('output', help='Output text file')
    
    args = parser.parse_args()
    
    try:
        extract_gff_info(args.input, args.output)
        print(f"Successfully extracted data from {args.input} to {args.output}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
        
    return 0

if __name__ == "__main__":
    sys.exit(main())
