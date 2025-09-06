#!/bin/bash

# scRNA-seq Data Processing Pipeline Script
# This script automates:
# 1. Merging R1 and R2 files
# 2. Aligning using HISAT2
# 3. Converting SAM to BAM
# 4. Sorting BAM files
# 5. Cleaning up intermediate files
# 6. Parsing BAM file to generate unique reads and duplicate reads reports

# Exit on error
set -e

# Default parameters
THREADS=16
CHUNK_SIZE=20000
MAX_INTRON_LEN=5000
MIN_READS=5000  # Minimum number of reads per barcode
GENOME_INDEX="/home2/tvieira/genomes/Pvivax-Sal1/Pvivax-Sal1_index"
BAM_PARSER="/home2/tvieira/scripts/stream_bam_parsing.py"  # Path to bam_parsing.py script

# Function to display usage information
usage() {
    echo "Usage: $0 -1 <r1_file.fastq.gz> -2 <r2_file.fastq.gz> [options]"
    echo ""
    echo "Required arguments:"
    echo "  -1, --r1          R1 FASTQ file (gzipped)"
    echo "  -2, --r2          R2 FASTQ file (gzipped)"
    echo ""
    echo "Optional arguments:"
    echo "  -o, --output      Output prefix (default: derived from R1 filename)"
    echo "  -t, --threads     Number of threads (default: 8)"
    echo "  -c, --chunk-size  Chunk size for merging (default: 20000)"
    echo "  -g, --genome      Path to genome index (default: $GENOME_INDEX)"
    echo "  -m, --max-intron  Maximum intron length (default: 5000)"
    echo "  -p, --parser      Path to BAM parser script (default: $BAM_PARSER)"
    echo "  -n, --min-reads   Minimum reads per barcode (default: $MIN_READS)"
    echo "  -h, --help        Display this help message"
    exit 1
}

# Parse command line arguments
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -1|--r1)
            R1_FILE="$2"
            shift 2
            ;;
        -2|--r2)
            R2_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_PREFIX="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -c|--chunk-size)
            CHUNK_SIZE="$2"
            shift 2
            ;;
        -g|--genome)
            GENOME_INDEX="$2"
            shift 2
            ;;
        -m|--max-intron)
            MAX_INTRON_LEN="$2"
            shift 2
            ;;
        -p|--parser)
            BAM_PARSER="$2"
            shift 2
            ;;
        -n|--min-reads)
            MIN_READS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

# Check for required arguments
if [ -z "$R1_FILE" ] || [ -z "$R2_FILE" ]; then
    echo "Error: Both R1 and R2 files are required."
    usage
fi

# Set output prefix if not specified
if [ -z "$OUTPUT_PREFIX" ]; then
    # Extract base name from R1 file without extensions
    OUTPUT_PREFIX=$(basename "$R1_FILE" | sed 's/\.fastq\.gz$//' | sed 's/\.fq\.gz$//' | sed 's/_R1//' | sed 's/_1//')
fi

# Create log directory
mkdir -p logs

echo "======== scRNA-seq Processing Pipeline ========"
echo "Starting processing for: $OUTPUT_PREFIX"
echo "R1 file: $R1_FILE"
echo "R2 file: $R2_FILE"
echo "Using $THREADS threads"
echo "=============================================="

# Step 1: Merge R1 and R2
echo "[$(date)] Step 1: Merging R1 and R2 files..."
MERGED_FASTQ="${OUTPUT_PREFIX}_merged.fastq"
python -u /home2/tvieira/scripts/merge_r1_r2.py "$R1_FILE" "$R2_FILE" \
    --threads "$THREADS" \
    --chunk-size "$CHUNK_SIZE" \
    --output "$MERGED_FASTQ" \
    &> "logs/${OUTPUT_PREFIX}_merge.log"
echo "[$(date)] Merging completed."

# Step 2: Align using HISAT2
echo "[$(date)] Step 2: Aligning merged reads with HISAT2..."
SAM_FILE="${OUTPUT_PREFIX}_mapped.sam"
ALIGNMENT_LOG="logs/${OUTPUT_PREFIX}_alignment.log"
hisat2 -p "$THREADS" \
    -x "$GENOME_INDEX" \
    -U "$MERGED_FASTQ" \
    -S "$SAM_FILE" \
    --no-unal \
    --max-intronlen "$MAX_INTRON_LEN" \
    &> "$ALIGNMENT_LOG"

# Keep only the last 6 lines of the alignment log
if [ -f "$ALIGNMENT_LOG" ]; then
    # Create a temporary file
    TEMP_LOG=$(mktemp)
    # Extract the last 6 lines
    tail -n 6 "$ALIGNMENT_LOG" > "$TEMP_LOG"
    # Replace the original log with the truncated version
    mv "$TEMP_LOG" "$ALIGNMENT_LOG"
fi

echo "[$(date)] Alignment completed."

# Step 3: Convert SAM to BAM
echo "[$(date)] Step 3: Converting SAM to BAM..."
BAM_FILE="${OUTPUT_PREFIX}_mapped.bam"
samtools view -@ "$THREADS" -b "$SAM_FILE" > "$BAM_FILE"
echo "[$(date)] SAM to BAM conversion completed."

# Step 4: Sort BAM file
echo "[$(date)] Step 4: Sorting BAM file..."
SORTED_BAM="${OUTPUT_PREFIX}_sorted.bam"
samtools sort -@ "$THREADS" -o "$SORTED_BAM" "$BAM_FILE"
echo "[$(date)] BAM sorting completed."

# Generate index for the sorted BAM (optional but recommended)
#echo "[$(date)] Step 5: Indexing sorted BAM file..."
#samtools index "$SORTED_BAM"
#echo "[$(date)] Indexing completed."

# Optionally remove intermediate files to save space
echo "[$(date)] Step 5: Cleaning up intermediate files..."
#if [ -f "$SORTED_BAM" ] && [ -f "${SORTED_BAM}.bai" ]; then
if [ -f "$SORTED_BAM" ]; then
    # Only remove if final products exist
    rm -f "$SAM_FILE"
    rm -f "$BAM_FILE"
    # Uncomment to remove merged FASTQ if desired
    rm -f "$MERGED_FASTQ"
    echo "Intermediate files removed."
else
    echo "Warning: Final BAM missing, keeping intermediate files."
fi

## There's no need to run this if you have multiple lanes
## merge_lanes_and_parse will merge the sorted BAM, 
## parse it and create the duplicate reads files

# Step 6: Run BAM parsing to generate read statistics
#echo "[$(date)] Step 6: Parsing BAM file for read statistics..."
#if [ -f "$SORTED_BAM" ] && [ -f "$BAM_PARSER" ]; then
#    PARSE_LOG="logs/${OUTPUT_PREFIX}_bam_parsing.log"
#    echo "Running BAM parser with $THREADS threads..."
#    
#    # Extract fields from sorted BAM and pipe to the parser script
#    samtools view -@ "$THREADS" "$SORTED_BAM" | \
#    python "$BAM_PARSER" \
#        --threads "$THREADS" \
#        --output "$OUTPUT_PREFIX" \
#        --chunk-size 100000 \
#        &> "$PARSE_LOG"
#    
#    echo "[$(date)] BAM parsing completed."
#    
#    # Check if output files were created
#    if [ -f "${OUTPUT_PREFIX}_unique_reads_per_cell.txt" ] && [ -f "${OUTPUT_PREFIX}_duplicate_reads.txt" ]; then
#        echo "BAM parsing successful - output files created."
#    else
#        echo "Warning: BAM parsing may have failed - output files not found."
#    fi
#else
#    echo "Warning: Cannot run BAM parsing. Either sorted BAM file or parser script is missing."
#    echo "BAM file: $SORTED_BAM"
#    echo "Parser script: $BAM_PARSER"
#fi
#
# Step 7: Filter barcodes based on minimum read count
#echo "[$(date)] Step 7: Filtering barcodes with fewer than $MIN_READS reads..."
#if [ -f "${OUTPUT_PREFIX}_unique_reads_per_cell.txt" ]; then
#    FILTERED_OUTPUT="${OUTPUT_PREFIX}_filtered_barcodes.txt"
#    FILTER_LOG="logs/${OUTPUT_PREFIX}_filtering.log"
#    
#    # Use awk to filter out barcodes with fewer than MIN_READS reads
#    # Skip the header line (NR > 1) and print barcodes with read count >= MIN_READS
#    awk -v min="$MIN_READS" 'NR == 1 {print $0; next} $2 >= min {print $0}' "${OUTPUT_PREFIX}_unique_reads_per_cell.txt" > "$FILTERED_OUTPUT" 2> "$FILTER_LOG"
    
    # Count the number of filtered barcodes
#    TOTAL_BARCODES=$(wc -l < "${OUTPUT_PREFIX}_unique_reads_per_cell.txt")
#    FILTERED_BARCODES=$(wc -l < "$FILTERED_OUTPUT")
#    TOTAL_BARCODES=$((TOTAL_BARCODES - 1))  # Subtract header line
#    FILTERED_BARCODES=$((FILTERED_BARCODES - 1))  # Subtract header line
    
#    echo "[$(date)] Barcode filtering completed."
#    echo "Original number of barcodes: $TOTAL_BARCODES"
#    echo "Barcodes after filtering (>= $MIN_READS reads): $FILTERED_BARCODES"
#    echo "Removed barcodes: $((TOTAL_BARCODES - FILTERED_BARCODES))"
#else
#    echo "Warning: Cannot filter barcodes. Unique reads file not found."
#fi

# Print summary
echo "============= Pipeline Complete =============="
echo "Output files:"
echo "  Sorted BAM: $SORTED_BAM"
if [ -f "${OUTPUT_PREFIX}_unique_reads_per_cell.txt" ]; then
    echo "  Unique reads per cell: ${OUTPUT_PREFIX}_unique_reads_per_cell.txt"
fi
if [ -f "${OUTPUT_PREFIX}_duplicate_reads.txt" ]; then
    echo "  Duplicate reads: ${OUTPUT_PREFIX}_duplicate_reads.txt"
fi
if [ -f "${OUTPUT_PREFIX}_filtered_barcodes.txt" ]; then
    echo "  Filtered barcodes (>= $MIN_READS reads): ${OUTPUT_PREFIX}_filtered_barcodes.txt"
fi
echo "Logs are available in the logs directory."
echo "=============================================="
