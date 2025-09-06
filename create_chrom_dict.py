#!/usr/bin/env python3
import argparse
import sys
import json
import pickle
from collections import defaultdict

def create_chromosome_dict(annotation_file, output_file=None, format='pickle'):
    """
    Create a dictionary with chromosomes as keys and lists of (gene, position) tuples as values.
    
    Args:
        annotation_file: Path to the annotation file
        output_file: Path to save the dictionary (None = don't save)
        format: Output format ('pickle' or 'json')
    
    Returns:
        The chromosome dictionary
    """
    # Use defaultdict to avoid checking if chromosome exists before appending
    chr_dic = defaultdict(list)
    
    try:
        with open(annotation_file, "r") as ann:
            # First pass: gather all chromosome names
            chromosomes = set()
            for line in ann:
                parts = line.strip().split()
                if len(parts) < 4:  # Skip malformed lines
                    continue
                chromosomes.add(parts[0])
            
            print(f"Found {len(chromosomes)} chromosomes in annotation file")
            
        # Second pass: populate the dictionary with gene data
        with open(annotation_file, "r") as annotation_file:
            for line_num, line in enumerate(annotation_file, 1):
                try:
                    parts = line.strip().split()
                    if len(parts) < 5:  # Expect at least 5 columns
                        print(f"Warning: Line {line_num} has incomplete data, skipping")
                        continue
                        
                    ann_chromossome = parts[0]
                    ann_gene = parts[1]
                    
                    # Convert positions to integers with error handling
                    try:
                        start_pos = int(parts[2])
                        end_pos = int(parts[3])
                    except ValueError:
                        print(f"Warning: Line {line_num} has non-integer positions, skipping")
                        continue
                        
                    ann_position = (start_pos, end_pos)
                    new_gene = (ann_gene, ann_position)
                    chr_dic[ann_chromossome].append(new_gene)
                    
                except Exception as e:
                    print(f"Warning: Error processing line {line_num}: {e}")
                    continue
        
        # Convert defaultdict to regular dict for output
        chr_dic = dict(chr_dic)
        
        # Print summary
        for chr_name, genes in chr_dic.items():
            print(f"Chromosome {chr_name}: {len(genes)} genes")
        
        # Save dictionary if output file is specified
        if output_file:
            if format.lower() == 'pickle':
                with open(output_file, 'wb') as f:
                    pickle.dump(chr_dic, f)
                print(f"Dictionary saved as pickle to {output_file}")
            elif format.lower() == 'json':
                # Convert tuple positions to lists for JSON serialization
                json_dic = {}
                for chr_name, genes in chr_dic.items():
                    json_dic[chr_name] = []
                    for gene, pos in genes:
                        json_dic[chr_name].append([gene, list(pos)])
                
                with open(output_file, 'w') as f:
                    json.dump(json_dic, f, indent=2)
                print(f"Dictionary saved as JSON to {output_file}")
            else:
                print(f"Warning: Unsupported format '{format}', dictionary not saved")
        
        return chr_dic
        
    except FileNotFoundError:
        print(f"Error: Annotation file '{annotation_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Create chromosome dictionary from annotation file')
    parser.add_argument('input', help='Input annotation file')
    parser.add_argument('--output', '-o', help='Output file to save dictionary')
    parser.add_argument('--format', '-f', choices=['pickle', 'json'], default='pickle',
                        help='Output format (default: pickle)')
    
    args = parser.parse_args()
    
    # Create the dictionary
    chr_dic = create_chromosome_dict(args.input, args.output, args.format)
    
    print(f"Successfully created chromosome dictionary with {len(chr_dic)} chromosomes")
    return 0

if __name__ == "__main__":
    sys.exit(main())
