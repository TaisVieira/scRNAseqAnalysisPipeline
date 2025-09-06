import re
import threading
from queue import Queue
import gzip
import argparse
from concurrent.futures import ThreadPoolExecutor
import os

class FastqFormatError(Exception):
    """Custom exception for FASTQ format violations"""
    pass

class FastqChunkProcessor:
    def __init__(self, chunk_size=10000):
        self.chunk_size = chunk_size
        self.trimmed_counter = 0
        self.total_counter = 0
        self.lock = threading.Lock()

    def validate_chunk(self, chunk):
        """Validate that chunk contains proper FASTQ format"""
        if not chunk:
            return

        # Check if first line starts with @
        if not chunk[0].startswith('@'):
            raise FastqFormatError(f"Chunk doesn't start with '@': {chunk[0][:50]}...")

        # Check if chunk length is multiple of 4
        if len(chunk) % 4 != 0:
            raise FastqFormatError(f"Chunk length {len(chunk)} is not a multiple of 4")

        # Validate each read in chunk
        for i in range(0, len(chunk), 4):
            if not chunk[i].startswith('@'):
                raise FastqFormatError(f"Read header doesn't start with '@' at line {i}: {chunk[i][:50]}...")
            if not chunk[i+2].strip() == '+':
                raise FastqFormatError(f"Third line is not '+' at line {i+2}: {chunk[i+2][:50]}...")

        return True

    def process_chunk(self, r1_chunk, r2_chunk, output_queue):
        try:
            # Validate both chunks before processing
            self.validate_chunk(r1_chunk)
            self.validate_chunk(r2_chunk)

            local_output = []
            local_trimmed = 0
            local_total = 0

            chunk_lines = zip(r1_chunk, r2_chunk)
            line_counter = 0
            current_record = []
            i = 0

            for r1_line, r2_line in chunk_lines:
                if re.search('^@', r1_line) and re.search('^@', r2_line):
                    line_counter = 1
                    r1_id = list(r1_line.split(" "))[0]
                    current_id = [r1_id]
                else:
                    line_counter += 1
                    if line_counter == 2:
                        barcode = r1_line[:12]
                        umi = r1_line[12:]
                        new_header = current_id[0]+'_'+barcode+'_'+umi
                        mrna = r2_line

                        # Process polyA trimming
                        trim_counter = 0
                        for index, letter in enumerate(mrna):
                            if letter == "A":
                                trim_counter += 1
                                if trim_counter == 6:
                                    trim_counter = 0
                                    mrna = mrna[:index-6] + "\n"
                                    i = index
                                    local_trimmed += 1
                                    break
                            else:
                                trim_counter = 0

                        current_record.extend([new_header, mrna, "+\n"])

                    elif line_counter == 4:
                        if i == 0:
                            quality = r2_line
                        else:
                            quality = r2_line[:i-6] + "\n"
                            i = 0

                        current_record.append(quality)
                        local_output.extend(current_record)
                        current_record = []
                        line_counter = 0
                        local_total += 1

            output_queue.put((local_output, local_trimmed, local_total))
        except FastqFormatError as e:
            print(f"FASTQ format error in chunk: {str(e)}")
            raise
        except Exception as e:
            print(f"Error processing chunk: {str(e)}")
            raise

def read_in_chunks(file, chunk_size):
    chunk = []
    line_count = 0
    for line in file:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        chunk.append(line)
        line_count += 1
        if line_count == chunk_size * 4:  # Each record has 4 lines
            if not chunk[0].startswith('@'):
                raise FastqFormatError(f"Chunk doesn't start with '@': {chunk[0][:50]}...")
            yield chunk
            chunk = []
            line_count = 0
    if chunk:  # Don't forget the last chunk
        if not len(chunk) % 4 == 0:
            raise FastqFormatError(f"Final chunk has {len(chunk)} lines, not a multiple of 4")
        if not chunk[0].startswith('@'):
            raise FastqFormatError(f"Final chunk doesn't start with '@': {chunk[0][:50]}...")
        yield chunk


def main(r1_path, r2_path, output_path="output.fastq", num_threads=4, chunk_size=10000):
    import time
    import datetime

    # Initialize counter and start time
    chunk_counter = 0
    start_time = time.time()

    processor = FastqChunkProcessor(chunk_size)
    output_queue = Queue()

    print(f"Starting processing with {num_threads} threads...", flush=True)
    print(f"Processing files:", flush=True)
    print(f"R1: {r1_path}", flush=True)
    print(f"R2: {r2_path}", flush=True)
    print(f"Output: {output_path}", flush=True)

    # Get the directory for the output file
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    control_file_path = os.path.join(os.path.dirname(output_path), "read_merge_controls.txt")

    with open(output_path, "w") as output_file, open(control_file_path, "w") as control_file:
        opener = gzip.open if r1_path.endswith('.gz') else open
        with opener(r1_path, 'rt') as fp1, opener(r2_path, 'rt') as fp2:
            with ThreadPoolExecutor(max_workers=num_threads) as executor:
                futures = []
                try:
                    for r1_chunk, r2_chunk in zip(read_in_chunks(fp1, chunk_size),
                                                read_in_chunks(fp2, chunk_size)):
                        future = executor.submit(processor.process_chunk,
                                              r1_chunk, r2_chunk, output_queue)
                        futures.append(future)
                        chunk_counter += 1
                        if chunk_counter % 1000 == 0:
                            elapsed = time.time() - start_time
                            elapsed_minutes = elapsed / 60
                            print(f"[{datetime.datetime.now()}] Processed {chunk_counter} chunks. Time elapsed: {elapsed_minutes:.2f}m", flush=True)

                    # Wait for all threads to complete
                    for future in futures:
                        future.result()

                    # Process results
                    total_trimmed = 0
                    total_sequences = 0

                    print(f"[{datetime.datetime.now()}] All chunks processed. Writing output...", flush=True)

                    while not output_queue.empty():
                        chunk_output, chunk_trimmed, chunk_total = output_queue.get()
                        for line in chunk_output:
                            print(line, end='', file=output_file)
                        total_trimmed += chunk_trimmed
                        total_sequences += chunk_total


                    percent_trim = (total_trimmed / total_sequences) * 100
                    control = f"Total number of sequences: {total_sequences}\nTrimmed sequences: {total_trimmed}\nPercentage of trimmed sequences: {percent_trim:.2f}"
                    print(control, file=control_file)

                    elapsed = time.time() - start_time
                    elapsed_minutes = elapsed / 60
                    print(f"[{datetime.datetime.now()}] Processing completed successfully! Total time: {elapsed_minutes:.2f}m", flush=True)

                except FastqFormatError as e:
                    print(f"\n[{datetime.datetime.now()}] Error: Invalid FASTQ format detected - {str(e)}", flush=True)
                    raise
                except Exception as e:
                    print(f"\n[{datetime.datetime.now()}] Error during processing: {str(e)}", flush=True)
                    raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process FASTQ files with multiple threads')
    parser.add_argument('r1', help='Path to R1 FASTQ file')
    parser.add_argument('r2', help='Path to R2 FASTQ file')
    parser.add_argument('--output', help='Path to output FASTQ file', default='output.fastq')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads to use')
    parser.add_argument('--chunk-size', type=int, default=10000, help='Number of records per chunk')

    args = parser.parse_args()

    main(args.r1, args.r2, args.output, args.threads, args.chunk_size)
