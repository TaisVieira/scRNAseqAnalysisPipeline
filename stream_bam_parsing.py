#!/usr/bin/env python3
import argparse
import sys
import threading
import queue
import time
import os
import tempfile
import pickle
from collections import defaultdict

def process_chunk(chunk_lines, thread_id):
    """Process a chunk of input lines and return results."""
    # Only track duplicate counts since unique_reads_per_cell is no longer needed
    local_dup_dic = defaultdict(int)

    for i, line in enumerate(chunk_lines):
        #if i % 10000000 == 0 and i > 0:
        #    print(f"Thread {thread_id}: Processed {i} reads")

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
                # Count duplicates with strand info
                unique_str = f"{barcode}, {umi}, {chromossome}, {position}, Forward"
                local_dup_dic[unique_str] += 1
            elif strand == "16":
                strand = "Reverse"
                # Count duplicates with strand info
                unique_str = f"{barcode}, {umi}, {chromossome}, {position}, Reverse"
                local_dup_dic[unique_str] += 1
            else:
                continue  # Skip secondary alignments

        except Exception as e:
            pass  # Skip problematic lines silently for performance

    return local_dup_dic

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

def worker_thread(thread_id, queue_obj, temp_dir):
    """Process chunks from the queue and write results to temporary files."""
    temp_file_count = 0
    processed_chunks = 0

    while True:
        chunk = queue_obj.get()

        if chunk is None:  # End signal
            queue_obj.task_done()
            break

        local_dup_dic = process_chunk(chunk, thread_id)

        # Save results to temporary files
        temp_file_path = os.path.join(temp_dir, f"thread_{thread_id}_part_{temp_file_count}.pkl")
        with open(temp_file_path, 'wb') as f:
            pickle.dump(local_dup_dic, f)

        # Clear local dictionary to free memory
        local_dup_dic.clear()

        temp_file_count += 1
        processed_chunks += 1

        #if processed_chunks % 10 == 0:
        #   print(f"Thread {thread_id}: Saved {processed_chunks} chunks to temporary files")

        queue_obj.task_done()

    #print(f"Thread {thread_id}: Completed, saved {temp_file_count} temporary files")
    return temp_file_count

def stream_merge_results(temp_dir, output_prefix):
    """Stream through temporary files and directly write output in a memory-efficient way."""
    print("Starting streaming merge of results...")

    # For duplicate counts
    dup_counts = defaultdict(int)

    # Get all temporary files
    temp_files = [f for f in os.listdir(temp_dir) if f.endswith('.pkl')]
    total_files = len(temp_files)

    # Process temporary files one by one
    for i, temp_file in enumerate(temp_files):
        if (i + 1) % 300 == 0 or (i + 1) == total_files:
            print(f"Merging file {i + 1}/{total_files}")

        full_path = os.path.join(temp_dir, temp_file)

        try:
            with open(full_path, 'rb') as f:
                local_dup_dic = pickle.load(f)

                # Process duplicate dictionary
                for unique_str, count in local_dup_dic.items():
                    dup_counts[unique_str] += count

                # Clear loaded dictionaries to free memory
                local_dup_dic.clear()

            # Delete the temporary file after processing
            os.remove(full_path)

        except Exception as e:
            print(f"Error processing temporary file {temp_file}: {e}")

    # Write output file directly from accumulated counts
    write_streaming_output(dup_counts, output_prefix)

    # Clear dictionaries
    dup_counts.clear()

    print("Streaming merge completed")

def write_streaming_output(dup_counts, output_prefix):
    """Write results to output file from the streaming process."""
    print(f"Writing output to {output_prefix}_duplicate_reads.txt")

    # Only write the duplicate reads file since unique_reads_per_cell is no longer needed
    with open(f"{output_prefix}_duplicate_reads.txt", "w+") as duplicate:
        print('Mapped_sequences\tDuplicates', file=duplicate)
        for unique_str, count in dup_counts.items():
            print(f"{unique_str}\t{count}", file=duplicate)

    # Report total reads processed
    total_reads = sum(dup_counts.values())
    print(f"Total reads processed: {total_reads}")

def main():
    parser = argparse.ArgumentParser(description='Process BAM data from stdin with multiple threads')
    parser.add_argument('--threads', '-t', type=int, default=4, help='Number of threads (default: 4)')
    parser.add_argument('--output', '-o', default='output', help='Output file prefix (default: output)')
    parser.add_argument('--chunk-size', '-c', type=int, default=100000, help='Chunk size (default: 100000 lines)')
    parser.add_argument('--temp-dir', '-d', default=None, help='Temp directory for intermediate files (default: system temp)')
    args = parser.parse_args()

    # Create temporary directory if not specified
    temp_dir = args.temp_dir
    if temp_dir is None:
        temp_dir = tempfile.mkdtemp(prefix='bam_processor_')
    else:
        os.makedirs(temp_dir, exist_ok=True)

    print(f"Using temporary directory: {temp_dir}")

    start_time = time.time()
    print(f"Starting processing with {args.threads} threads")

    # Create thread-safe queue
    chunk_queue = queue.Queue(maxsize=100)  # Limit queue size to avoid memory issues

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
            args=(i, chunk_queue, temp_dir)
        )
        worker.daemon = True
        worker.start()
        workers.append(worker)

    # Wait for all work to complete
    reader.join()
    for worker in workers:
        worker.join()

    # Stream merge results and write output
    stream_merge_results(temp_dir, args.output)

    # Clean up temporary directory if it was auto-created
    if args.temp_dir is None:
        try:
            os.rmdir(temp_dir)  # Only removes if empty
            print(f"Removed temporary directory: {temp_dir}")
        except OSError:
            print(f"Note: Could not remove temporary directory (it may not be empty): {temp_dir}")

    elapsed_time = time.time() - start_time
    elapsed_time_minutes = elapsed_time / 60
    print(f"Processing completed in {elapsed_time_minutes:.2f} minutes")

    return 0

if __name__ == "__main__":
    exit(main())
