#!/usr/bin/env python3
import pickle
import sys
import argparse
import csv
from collections import defaultdict, Counter

def map_reads_to_genes(duplicate_reads_file, chr_dict_file, min_reads=5000, output_file=None):
    """
    Map reads to genes using a chromosome dictionary and create a hit dictionary.
    
    Args:
        duplicate_reads_file: Path to the duplicate reads file
        chr_dict_file: Path to the pickled chromosome dictionary
        min_reads: Minimum number of reads for a cell to be considered
        output_file: Path to save the results in Seurat-compatible format
    
    Returns:
        The gene-cell counts dictionary
    """
    # Load the chromosome dictionary from pickle file
    try:
        with open(chr_dict_file, 'rb') as f:
            chr_dic = pickle.load(f)
        print(f"Loaded chromosome dictionary with {len(chr_dic)} chromosomes")
    except Exception as e:
        print(f"Error loading chromosome dictionary: {e}")
        sys.exit(1)
    
    # First pass to count reads per cell (if needed)
    cell_reads = {}
    if min_reads > 0:
        print(f"Counting reads per cell (minimum threshold: {min_reads})...")
        try:
            with open(duplicate_reads_file, "r") as mapped:
                next(mapped)  # Skip header
                for line_num, line in enumerate(mapped, 1):
                    try:
                        parts = line.strip().split('\t')
                        if len(parts) < 2:
                            continue
                        
                        mapped_sequence = parts[0].strip()
                        seq_parts = mapped_sequence.split(', ')
                        if len(seq_parts) < 5:
                            continue
                        
                        read_bc = seq_parts[0]
                        duplicates = int(parts[1])
                        
                        # Count the reads (including duplicates)
                        cell_reads[read_bc] = cell_reads.get(read_bc, 0) + duplicates
                        
                        if line_num % 100000000 == 0:
                            print(f"Processed {line_num:,} lines for cell counting")
                    except Exception as e:
                        if line_num % 100000000 == 0:
                            print(f"Warning at line {line_num}: {e}")
                        continue
            
            # Filter cells with fewer than min_reads
            qualified_cells = {cell: count for cell, count in cell_reads.items() if count >= min_reads}
            print(f"Found {len(qualified_cells)} cells with at least {min_reads} reads out of {len(cell_reads)} total cells")
        except Exception as e:
            print(f"Error during read counting: {e}")
            sys.exit(1)
    else:
        # If min_reads is 0 or negative, consider all cells
        qualified_cells = None
        print("Skipping read count filtering (all cells will be considered)")
    
    # Create the hit dictionary and count UMIs per gene-cell combination
    gene_cell_counts = defaultdict(Counter)
    processed_reads = 0
    mapped_reads = 0
    
    try:
        with open(duplicate_reads_file, "r") as mapped:
            next(mapped)  # Skip header
            for line_num, line in enumerate(mapped, 1):
                try:
                    parts = line.strip().split('\t')
                    if len(parts) < 2:
                        continue
                    
                    mapped_sequence = parts[0].strip()
                    seq_parts = mapped_sequence.split(', ')
                    if len(seq_parts) < 5:
                        continue
                    
                    read_bc = seq_parts[0]
                    umi = seq_parts[1]  # Not used directly but available
                    read_chromossome = seq_parts[2]
                    
                    try:
                        read_position = int(seq_parts[3])
                        duplicates = int(parts[1])
                    except ValueError:
                        continue
                    
                    # Skip if the cell doesn't have enough reads
                    if qualified_cells is not None and read_bc not in qualified_cells:
                        continue
                    
                    processed_reads += duplicates
                    if line_num % 10000000 == 0:
                        print(f"Processed {line_num:,} lines ({processed_reads:,} reads), mapped {mapped_reads:,} to genes")
                    
                    # Look for the read in the chromosome dictionary
                    if read_chromossome in chr_dic:
                        found_match = False
                        for name, position in chr_dic[read_chromossome]:
                            pos1, pos2 = position
                            if pos1 < read_position < pos2:
                                # Increment counter by duplicate count
                                gene_cell_counts[name][read_bc] += duplicates
                                mapped_reads += duplicates
                                found_match = True
                                break  # No need to check other genes once we found a match
                        
                        if not found_match and line_num % 10000000 == 0:
                            print(f"No gene match for position {read_position} on {read_chromossome}")
                    
                except Exception as e:
                    if line_num % 1000000 == 0:  # Only show occasional warnings
                        print(f"Warning at line {line_num}: {e}")
                    continue
    
    except Exception as e:
        print(f"Error during read mapping: {e}")
        sys.exit(1)
    
    print(f"Finished processing {processed_reads:,} reads from {line_num:,} lines")
    print(f"Mapped {mapped_reads:,} reads to {len(gene_cell_counts)} genes")
    
    # Get all unique cell barcodes from the mapping
    all_cells = set()
    for gene_counts in gene_cell_counts.values():
        all_cells.update(gene_counts.keys())
    
    print(f"Found {len(all_cells)} unique cells in the mapped data")
    
    # Save results in Seurat-compatible format
    if output_file:
        try:
            # Create a matrix file (genes x cells)
            with open(output_file, 'w', newline='') as f:
                # Prepare header: gene name in first column, followed by cell barcodes
                cell_barcodes = sorted(list(all_cells))
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(['Gene'] + cell_barcodes)
                
                # Write each gene's expression counts across all cells
                for gene, cell_counts in sorted(gene_cell_counts.items()):
                    row = [gene]
                    for cell in cell_barcodes:
                        # Add the count for this gene-cell combination (0 if no reads mapped)
                        row.append(cell_counts.get(cell, 0))
                    writer.writerow(row)
            
            print(f"Saved Seurat-compatible matrix to {output_file}")
            
            # Also save a pickle version of the gene-cell counts for future use
            pickle_file = output_file + ".pkl"
            with open(pickle_file, 'wb') as f:
                pickle.dump(dict(gene_cell_counts), f)
            print(f"Also saved raw counts data to {pickle_file}")
            
        except Exception as e:
            print(f"Error saving output: {e}")
    
    return gene_cell_counts

def main():
    parser = argparse.ArgumentParser(description='Map reads to genes and prepare data for Seurat')
    parser.add_argument('duplicate_reads', help='Duplicate reads file')
    parser.add_argument('chr_dict', help='Pickled chromosome dictionary file')
    parser.add_argument('--min-reads', '-m', type=int, default=5000,
                       help='Minimum reads per cell (default: 5000, use 0 to include all cells)')
    parser.add_argument('--output', '-o', required=True, help='Output file for Seurat-compatible matrix')
    
    args = parser.parse_args()
    
    # Map reads to genes
    gene_cell_counts = map_reads_to_genes(args.duplicate_reads, args.chr_dict, args.min_reads, args.output)
    
    # Print some statistics about the results
    total_umi_counts = sum(sum(cell_counts.values()) for cell_counts in gene_cell_counts.values())
    print(f"\nTotal UMI counts in matrix: {total_umi_counts:,}")
    
    # Find genes with highest expression
    gene_total_counts = {gene: sum(cell_counts.values()) for gene, cell_counts in gene_cell_counts.items()}
    
    print("\nTop 5 genes by total UMI counts:")
    for gene, count in sorted(gene_total_counts.items(), key=lambda x: x[1], reverse=True)[:5]:
        print(f"{gene}: {count:,} UMIs")
    
    # Find cells with most genes detected
    #cell_gene_counts = defaultdict(int)
    #for gene, cell_counts in gene_cell_counts.items():
    #    for cell, count in cell_counts.items():
    #        if count > 0:
    #            cell_gene_counts[cell] += 1
    
    #if cell_gene_counts:
    #    print("\nTop 5 cells by genes detected:")
    #    for cell, count in sorted(cell_gene_counts.items(), key=lambda x: x[1], reverse=True)[:5]:
    #        print(f"{cell}: {count} genes")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
