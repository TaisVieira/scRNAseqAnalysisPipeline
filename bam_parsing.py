#!/usr/bin/env python3
import argparse
import sys
import threading
import queue
import time
from collections import defaultdict

def process_chunk(chunk_lines, thread_id):
    """Process a chunk of input lines and return results."""
    local_read_dic = defaultdict(list)
    local_dup_dic = defaultdict(int)
    
    for i, line in enumerate(chunk_lines):
        if i % 10000000 == 0 and i > 0:
            print(f"Thread {thread_id}: Processed {i} reads")
            
        try:
            info = line.split("\t")
            read = info[0]
            strand = info[1]
            
            # Skip if unable to parse
            if len(info) < 4 or not read or '_' not in read:
                continue
                
            # Extract information
            barcode = read.split("_")[1]
            umi = read.split("_")[2]
            chromossome = info[2]
            position = info[3]
            
            # Check if it's a primary alignment
            if strand == "0":
                strand = "Forward"
                unique = [barcode, umi, chromossome, position, strand]
            elif strand == "16":
                strand = "Reverse"
                unique = [barcode, umi, chromossome, position, strand]
            else:
                continue  # Skip secondary alignments
                
            # Store unique reads per cell
            if unique not in local_read_dic[barcode]:
                local_read_dic[barcode].append(unique)
                
            # Count duplicates
            unique_str = f"{barcode}, {umi}, {chromossome}, {position}"
            local_dup_dic[unique_str] += 1
            
        except Exception as e:
            pass  # Skip problematic lines silently for performance
    
    return local_read_dic, local_dup_dic

def reader_thread(input_stream, chunk_size, queue_obj, num_threads):
    """Read input and distribute chunks to worker threads."""
    total_lines = 0
    current_chunk = []
    
    try:
        for line in input_stream:
            current_chunk.append(line)
            total_lines += 1
            
            if len(current_chunk) >= chunk_size:
                queue_obj.put(current_chunk)
                current_chunk = []
                
            if total_lines % 20000000 == 0:
                print(f"Read {total_lines} lines so far")
                
    except Exception as e:
        print(f"Error in reader: {e}")
    
    # Put any remaining lines
    if current_chunk:
        queue_obj.put(current_chunk)
        
    # Add None signals to indicate end of data
    for _ in range(num_threads):
        queue_obj.put(None)
        
    print(f"Reader: Finished reading {total_lines} total lines")

def worker_thread(thread_id, queue_obj, results):
    """Process chunks from the queue."""
    while True:
        chunk = queue_obj.get()
        
        if chunk is None:  # End signal
            break
            
        result = process_chunk(chunk, thread_id)
        results.append(result)
        queue_obj.task_done()

def merge_results(results):
    """Merge results from all threads."""
    merged_read_dic = defaultdict(list)
    merged_dup_dic = defaultdict(int)
    
    print("Merging results from all threads...")
    for local_read_dic, local_dup_dic in results:
        for barcode, reads in local_read_dic.items():
            for read in reads:
                if read not in merged_read_dic[barcode]:
                    merged_read_dic[barcode].append(read)
        
        for unique_str, count in local_dup_dic.items():
            merged_dup_dic[unique_str] += count
    
    return merged_read_dic, merged_dup_dic

def write_output(read_dic, dup_dic, output_prefix):
    """Write results to output files."""
    print(f"Writing output to {output_prefix}_unique_reads_per_cell.txt and {output_prefix}_duplicate_reads.txt")
    
    with open(f"{output_prefix}_unique_reads_per_cell.txt", "w+") as unique_cell:
        print('Barcode\tReads', file=unique_cell)
        for k, v in read_dic.items():
            print(f"{k}\t{len(v)}", file=unique_cell)
                  
    with open(f"{output_prefix}_duplicate_reads.txt", "w+") as duplicate:
        print('Mapped_sequences\tDuplicates', file=duplicate)
        for k, v in dup_dic.items():
            print(f"{k}\t{v}", file=duplicate)

def main():
    parser = argparse.ArgumentParser(description='Process BAM data from stdin with multiple threads')
    parser.add_argument('--threads', '-t', type=int, default=4, help='Number of threads (default: 4)')
    parser.add_argument('--output', '-o', default='output', help='Output file prefix (default: output)')
    parser.add_argument('--chunk-size', '-c', type=int, default=100000, help='Chunk size (default: 100000 lines)')
    args = parser.parse_args()
    
    start_time = time.time()
    print(f"Starting processing with {args.threads} threads")
    
    # Create thread-safe queue and results list
    chunk_queue = queue.Queue(maxsize=100)  # Limit queue size to avoid memory issues
    results = []
    
    # Start reader thread to read from stdin
    reader = threading.Thread(
        target=reader_thread, 
        args=(sys.stdin, args.chunk_size, chunk_queue, args.threads)
    )
    reader.daemon = True
    reader.start()
    
    # Start worker threads
    workers = []
    for i in range(args.threads):
        worker = threading.Thread(
            target=worker_thread,
            args=(i, chunk_queue, results)
        )
        worker.daemon = True
        worker.start()
        workers.append(worker)
    
    # Wait for all work to complete
    reader.join()
    for worker in workers:
        worker.join()
    
    # Merge and write results
    read_dic, dup_dic = merge_results(results)
    write_output(read_dic, dup_dic, args.output)
    
    elapsed_time = time.time() - start_time
    print(f"Processing completed in {elapsed_time:.2f} seconds")
    print(f"Found {sum(len(v) for v in read_dic.values())} unique reads")
    print(f"Found {sum(v for v in dup_dic.values())} total reads")
    
    return 0

if __name__ == "__main__":
    exit(main())
