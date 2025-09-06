#!/bin/bash

# Merge BAM files from multiple lanes and run stream_bam_parsing
# This script takes lane directories as input, finds BAM files, merges them,
# and runs the stream BAM parsing pipeline

set -e

# Default parameters
THREADS=16
CHUNK_SIZE=100000
TEMP_DIR=""
OUTPUT_PREFIX=""
BAM_PARSER="/home2/tvieira/scripts/stream_bam_parsing.py"
CLEANUP=true
BAM_PATTERN="*sorted.bam"  # Pattern to match BAM files

# Function to display usage
usage() {
    echo "Usage: $0 [options] <lane_dir1> <lane_dir2> [lane_dir3] ..."
    echo ""
    echo "This script merges BAM files from multiple lane directories and runs stream BAM parsing."
    echo ""
    echo "Arguments:"
    echo "  lane_dir1, lane_dir2, ...  Directories containing BAM files to merge"
    echo ""
    echo "Options:"
    echo "  -o, --output       Output prefix for final results (default: merged_lanes)"
    echo "  -t, --threads      Number of threads (default: 16)"
    echo "  -c, --chunk-size   Chunk size for BAM parsing (default: 100000)"
    echo "  -d, --temp-dir     Temporary directory for intermediate files (default: auto)"
    echo "  -p, --parser       Path to stream_bam_parsing.py script (default: $BAM_PARSER)"
    echo "  -b, --bam-pattern  Pattern to match BAM files (default: $BAM_PATTERN)"
    echo "  --no-cleanup       Don't remove intermediate merged BAM file"
    echo "  -h, --help         Display this help message"
    echo ""
    echo "Example:"
    echo "  $0 -o sample1_merged -t 20 lane1/ lane2/ lane3/"
    exit 1
}

# Parse command line arguments
LANE_DIRS=()
while [[ $# -gt 0 ]]; do
    case $1 in
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
        -d|--temp-dir)
            TEMP_DIR="$2"
            shift 2
            ;;
        -p|--parser)
            BAM_PARSER="$2"
            shift 2
            ;;
        -b|--bam-pattern)
            BAM_PATTERN="$2"
            shift 2
            ;;
        --no-cleanup)
            CLEANUP=false
            shift
            ;;
        -h|--help)
            usage
            ;;
        -*)
            echo "Unknown option: $1"
            usage
            ;;
        *)
            LANE_DIRS+=("$1")
            shift
            ;;
    esac
done

# Check if at least one lane directory is provided
if [ ${#LANE_DIRS[@]} -eq 0 ]; then
    echo "Error: At least one lane directory must be provided."
    usage
fi

# Set default output prefix if not specified
if [ -z "$OUTPUT_PREFIX" ]; then
    OUTPUT_PREFIX="merged_lanes"
fi

# Set default temp directory if not specified
if [ -z "$TEMP_DIR" ]; then
    TEMP_DIR="temp_bam_parsing"
fi

# Create logs directory
mkdir -p logs

echo "======== Lane Merging and BAM Parsing Pipeline ========"
echo "Lane directories: ${LANE_DIRS[*]}"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Using $THREADS threads"
echo "BAM pattern: $BAM_PATTERN"
echo "Temp directory: $TEMP_DIR"
echo "======================================================"

# Step 1: Find all BAM files in the specified directories
echo "[$(date)] Step 1: Finding BAM files in lane directories..."
BAM_FILES=()
for lane_dir in "${LANE_DIRS[@]}"; do
    if [ ! -d "$lane_dir" ]; then
        echo "Warning: Directory $lane_dir does not exist, skipping..."
        continue
    fi
    
    echo "Searching for BAM files in $lane_dir..."
    
    # Find BAM files matching the pattern
    while IFS= read -r -d '' bam_file; do
        if [ -f "$bam_file" ]; then
            BAM_FILES+=("$bam_file")
            echo "  Found: $bam_file"
        fi
    done < <(find "$lane_dir" -name "$BAM_PATTERN" -type f -print0)
done

# Check if any BAM files were found
if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "Error: No BAM files found matching pattern '$BAM_PATTERN' in the specified directories."
    echo "Please check your lane directories and BAM file pattern."
    exit 1
fi

echo "Found ${#BAM_FILES[@]} BAM files total:"
printf '  %s\n' "${BAM_FILES[@]}"

# Step 2: Merge BAM files
echo "[$(date)] Step 2: Merging BAM files..."
MERGED_BAM="${OUTPUT_PREFIX}_merged.bam"

if [ ${#BAM_FILES[@]} -eq 1 ]; then
    echo "Only one BAM file found, creating symlink instead of merging..."
    ln -sf "$(realpath "${BAM_FILES[0]}")" "$MERGED_BAM"
else
    echo "Merging ${#BAM_FILES[@]} BAM files into $MERGED_BAM..."
    
    # Use samtools merge to combine all BAM files
    # The -@ option specifies the number of threads
    samtools merge -@ "$THREADS" "$MERGED_BAM" "${BAM_FILES[@]}" \
        > "logs/${OUTPUT_PREFIX}_merge.log" 2>&1
    
    echo "BAM files merged successfully."
fi

# Verify merged BAM file exists and is not empty
if [ ! -f "$MERGED_BAM" ] || [ ! -s "$MERGED_BAM" ]; then
    echo "Error: Merged BAM file is missing or empty: $MERGED_BAM"
    exit 1
fi

# Step 3: Run stream BAM parsing
echo "[$(date)] Step 3: Running stream BAM parsing..."

if [ ! -f "$BAM_PARSER" ]; then
    echo "Error: BAM parser script not found: $BAM_PARSER"
    exit 1
fi

# Create temp directory for BAM parsing
mkdir -p "$TEMP_DIR"

echo "Starting BAM parsing with the following parameters:"
echo "  Threads: $THREADS"
echo "  Chunk size: $CHUNK_SIZE"
echo "  Temp directory: $TEMP_DIR"
echo "  Output prefix: $OUTPUT_PREFIX"

# Run the BAM parsing pipeline
PARSE_LOG="logs/${OUTPUT_PREFIX}_bam_parsing.log"
echo "Running: samtools view -@ $THREADS $MERGED_BAM | python -u $BAM_PARSER -t $THREADS -o $OUTPUT_PREFIX -c $CHUNK_SIZE -d $TEMP_DIR"

samtools view -@ "$THREADS" "$MERGED_BAM" | \
    python -u "$BAM_PARSER" \
        -t "$THREADS" \
        -o "$OUTPUT_PREFIX" \
        -c "$CHUNK_SIZE" \
        -d "$TEMP_DIR" \
        > "$PARSE_LOG" 2>&1

echo "[$(date)] BAM parsing completed."

# Step 4: Cleanup intermediate files
if [ "$CLEANUP" = true ]; then
    echo "[$(date)] Step 4: Cleaning up intermediate files..."
    
    # Only remove temp directory if empty - KEEP the merged BAM file
    if [ -d "$TEMP_DIR" ]; then
        rmdir "$TEMP_DIR" 2>/dev/null || echo "Note: Temp directory not empty, keeping: $TEMP_DIR"
    fi
    
    echo "Cleanup completed (merged BAM file preserved)."
else
    echo "[$(date)] Step 4: Skipping cleanup (--no-cleanup specified)"
    echo "All files preserved:"
    echo "  Merged BAM: $MERGED_BAM"
    echo "  Temp directory: $TEMP_DIR"
fi

# Step 5: Verify output files and print summary
echo "============= Pipeline Complete =============="

# Check if expected output files exist
if [ -f "${OUTPUT_PREFIX}_duplicate_reads.txt" ]; then
    echo "✓ Output file created: ${OUTPUT_PREFIX}_duplicate_reads.txt"
    
    # Count number of unique sequences and total reads
    UNIQUE_SEQUENCES=$(tail -n +2 "${OUTPUT_PREFIX}_duplicate_reads.txt" | wc -l)
    TOTAL_READS=$(tail -n +2 "${OUTPUT_PREFIX}_duplicate_reads.txt" | awk -F'\t' '{sum += $2} END {print sum}')
    
    echo "  └─ Unique mapped sequences: $UNIQUE_SEQUENCES"
    echo "  └─ Total reads processed: $TOTAL_READS"
else
    echo "✗ Warning: Expected output file not found: ${OUTPUT_PREFIX}_duplicate_reads.txt"
fi

if [ -f "$MERGED_BAM" ]; then
    echo "✓ Merged BAM file preserved: $MERGED_BAM"
    BAM_SIZE=$(du -h "$MERGED_BAM" | cut -f1)
    echo "  └─ File size: $BAM_SIZE"
else
    echo "✗ Warning: Merged BAM file missing: $MERGED_BAM"
fi

echo ""
echo "Logs available:"
echo "  Merge log: logs/${OUTPUT_PREFIX}_merge.log"
echo "  BAM parsing log: logs/${OUTPUT_PREFIX}_bam_parsing.log"
echo ""
echo "Input summary:"
echo "  Lane directories processed: ${#LANE_DIRS[@]}"
echo "  BAM files merged: ${#BAM_FILES[@]}"
echo "=============================================="
