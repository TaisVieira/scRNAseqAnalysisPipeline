#!/bin/bash

# Genome Preparation Script for scRNA-seq Pipeline
# This script automates:
# 1. Downloading FASTA and GFF files
# 2. Building HISAT2 genome index
# 3. Extracting gene information from GFF
# 4. Creating chromosome dictionary for gene mapping

# Exit on error
set -e

# Default parameters
THREADS=16
SCRIPTS_DIR="/home2/tvieira/scripts"  # Update this to your scripts directory
CLEANUP=true

# Function to display usage information
usage() {
    echo "Usage: $0 --fasta <fasta_url> --gff <gff_url> --output <genome_name> [options]"
    echo ""
    echo "Required arguments:"
    echo "  --fasta           URL to genome FASTA file"
    echo "  --gff             URL to genome GFF file"
    echo "  --output          Output name for genome (e.g., 'Pknowlesi-68')"
    echo ""
    echo "Optional arguments:"
    echo "  --threads         Number of threads for HISAT2 indexing (default: 16)"
    echo "  --scripts-dir     Directory containing your Python scripts (default: $SCRIPTS_DIR)"
    echo "  --no-cleanup      Keep intermediate files (default: cleanup enabled)"
    echo "  --output-dir      Output directory (default: current directory)"
    echo "  --help            Display this help message"
    echo ""
    echo "Example:"
    echo "  $0 --fasta https://example.com/genome.fa.gz --gff https://example.com/annotations.gff3.gz --output MyGenome"
    exit 1
}

# Function to check if required tools are available
check_dependencies() {
    echo "Checking dependencies..."
    
    local missing_tools=()
    
    if ! command -v hisat2-build >/dev/null 2>&1; then
        missing_tools+=("hisat2-build")
    fi
    
    if ! command -v wget >/dev/null 2>&1 && ! command -v curl >/dev/null 2>&1; then
        missing_tools+=("wget or curl")
    fi
    
    if ! command -v python >/dev/null 2>&1 && ! command -v python3 >/dev/null 2>&1; then
        missing_tools+=("python")
    fi
    
    if [ ${#missing_tools[@]} -ne 0 ]; then
        echo "Error: Missing required tools: ${missing_tools[*]}"
        echo "Please install the missing tools and try again."
        exit 1
    fi
    
    echo "All required tools are available."
}

# Function to download files
download_file() {
    local url=$1
    local output_file=$2
    
    echo "Downloading $url to $output_file..."
    
    if command -v wget >/dev/null 2>&1; then
        wget -O "$output_file" "$url"
    elif command -v curl >/dev/null 2>&1; then
        curl -L -o "$output_file" "$url"
    else
        echo "Error: Neither wget nor curl is available for downloading files."
        exit 1
    fi
    
    # Verify download
    if [ ! -f "$output_file" ] || [ ! -s "$output_file" ]; then
        echo "Error: Failed to download $url"
        exit 1
    fi
    
    echo "Successfully downloaded $output_file"
}

# Function to decompress files if needed
decompress_if_needed() {
    local file=$1
    local output_name=$2
    
    if [[ "$file" == *.gz ]]; then
        echo "Decompressing $file..."
        gunzip -c "$file" > "$output_name"
        if [ "$CLEANUP" = true ]; then
            rm "$file"
        fi
        echo "$output_name"
    else
        echo "$file"
    fi
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --fasta)
            FASTA_URL="$2"
            shift 2
            ;;
        --gff)
            GFF_URL="$2"
            shift 2
            ;;
        --output)
            GENOME_NAME="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --scripts-dir)
            SCRIPTS_DIR="$2"
            shift 2
            ;;
        --no-cleanup)
            CLEANUP=false
            shift
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check for required arguments
if [ -z "$FASTA_URL" ] || [ -z "$GFF_URL" ] || [ -z "$GENOME_NAME" ]; then
    echo "Error: Missing required arguments."
    usage
fi

# Set output directory - create genome folder
if [ -n "$OUTPUT_DIR" ]; then
    GENOME_DIR="$OUTPUT_DIR/$GENOME_NAME"
else
    GENOME_DIR="$GENOME_NAME"
fi

mkdir -p "$GENOME_DIR"
cd "$GENOME_DIR"

# Check dependencies
check_dependencies

# Check if Python scripts exist
PYTHON_CMD="python3"
if ! command -v python3 >/dev/null 2>&1; then
    PYTHON_CMD="python"
fi

GFF_EXTRACT_SCRIPT="$SCRIPTS_DIR/gff_extract.py"
CREATE_DICT_SCRIPT="$SCRIPTS_DIR/create_chrom_dict.py"

if [ ! -f "$GFF_EXTRACT_SCRIPT" ]; then
    echo "Error: GFF extract script not found at $GFF_EXTRACT_SCRIPT"
    exit 1
fi

if [ ! -f "$CREATE_DICT_SCRIPT" ]; then
    echo "Error: Create chromosome dictionary script not found at $CREATE_DICT_SCRIPT"
    exit 1
fi

echo "========== Genome Preparation Pipeline =========="
echo "Genome name: $GENOME_NAME"
echo "FASTA URL: $FASTA_URL"
echo "GFF URL: $GFF_URL"
echo "Using $THREADS threads"
echo "Scripts directory: $SCRIPTS_DIR"
echo "=================================================="

# Create output directory structure
mkdir -p "logs"

# Step 1: Download FASTA file
echo "[$(date)] Step 1: Downloading genome FASTA file..."
# Keep original filename from URL for FASTA
FASTA_ORIGINAL_NAME=$(basename "$FASTA_URL" | sed 's/\.gz$//')
FASTA_DOWNLOADED="$FASTA_ORIGINAL_NAME"
if [[ "$FASTA_URL" == *.gz ]]; then
    FASTA_DOWNLOADED="${FASTA_ORIGINAL_NAME}.gz"
fi

download_file "$FASTA_URL" "$FASTA_DOWNLOADED"

# Decompress FASTA if needed but keep original name
FASTA_FILE=$(decompress_if_needed "$FASTA_DOWNLOADED" "$FASTA_ORIGINAL_NAME")

# Step 2: Download GFF file
echo "[$(date)] Step 2: Downloading genome GFF file..."
# Keep original filename from URL for GFF
GFF_ORIGINAL_NAME=$(basename "$GFF_URL" | sed 's/\.gz$//')
GFF_DOWNLOADED="$GFF_ORIGINAL_NAME"
if [[ "$GFF_URL" == *.gz ]]; then
    GFF_DOWNLOADED="${GFF_ORIGINAL_NAME}.gz"
fi

download_file "$GFF_URL" "$GFF_DOWNLOADED"

# Decompress GFF if needed but keep original name
GFF_FILE=$(decompress_if_needed "$GFF_DOWNLOADED" "$GFF_ORIGINAL_NAME")

# Step 3: Build HISAT2 index
echo "[$(date)] Step 3: Building HISAT2 genome index..."
INDEX_PREFIX="${GENOME_NAME}_index"
INDEX_LOG="logs/${GENOME_NAME}_hisat2_build.log"

hisat2-build -p "$THREADS" "$FASTA_FILE" "$INDEX_PREFIX" &> "$INDEX_LOG"

# Verify index was created
if ls "${INDEX_PREFIX}".*.ht2 1> /dev/null 2>&1; then
    echo "[$(date)] HISAT2 index built successfully."
else
    echo "Error: HISAT2 index build failed. Check $INDEX_LOG for details."
    exit 1
fi

# Step 4: Extract gene information from GFF
echo "[$(date)] Step 4: Extracting gene information from GFF..."
SIMPLIFIED_GFF="${GENOME_NAME}_simpl_gff.txt"
GFF_EXTRACT_LOG="logs/${GENOME_NAME}_gff_extract.log"

$PYTHON_CMD "$GFF_EXTRACT_SCRIPT" "$GFF_FILE" "$SIMPLIFIED_GFF" &> "$GFF_EXTRACT_LOG"

if [ ! -f "$SIMPLIFIED_GFF" ] || [ ! -s "$SIMPLIFIED_GFF" ]; then
    echo "Error: GFF extraction failed. Check $GFF_EXTRACT_LOG for details."
    exit 1
fi

echo "[$(date)] Gene information extracted successfully."

# Step 5: Create chromosome dictionary
echo "[$(date)] Step 5: Creating chromosome dictionary..."
CHROM_DICT="${GENOME_NAME}_chrm_dict.pkl"
DICT_CREATE_LOG="logs/${GENOME_NAME}_dict_create.log"

$PYTHON_CMD "$CREATE_DICT_SCRIPT" "$SIMPLIFIED_GFF" --output "$CHROM_DICT" --format pickle &> "$DICT_CREATE_LOG"

if [ ! -f "$CHROM_DICT" ]; then
    echo "Error: Chromosome dictionary creation failed. Check $DICT_CREATE_LOG for details."
    exit 1
fi

echo "[$(date)] Chromosome dictionary created successfully."

# Step 6: Create configuration file for easy reference
echo "[$(date)] Step 6: Creating configuration file..."
CONFIG_FILE="${GENOME_NAME}_config.txt"

cat > "$CONFIG_FILE" << EOF
# Genome Configuration for $GENOME_NAME
# Generated on $(date)

# File locations
GENOME_FASTA=$PWD/$FASTA_FILE
GENOME_GFF=$PWD/$GFF_FILE
HISAT2_INDEX=$PWD/$INDEX_PREFIX
CHROMOSOME_DICT=$PWD/$CHROM_DICT
SIMPLIFIED_GFF=$PWD/$SIMPLIFIED_GFF

# For use in your scripts:
# Update your process_scrna.sh with:
# GENOME_INDEX="$PWD/$INDEX_PREFIX"

# For gene mapping step, use:
# CHROMOSOME_DICT="$PWD/$CHROM_DICT"
EOF

echo "Configuration saved to $CONFIG_FILE"

# Optional cleanup
if [ "$CLEANUP" = true ]; then
    echo "[$(date)] Step 7: Cleaning up intermediate files..."
    # Keep the final processed files but remove any temporary downloads
    # The decompressed files are already handled in decompress_if_needed function
    echo "Cleanup completed."
fi

# Print summary
echo "=============== Genome Preparation Complete ==============="
echo "Genome name: $GENOME_NAME"
echo ""
echo "Generated files:"
echo "  HISAT2 Index: $INDEX_PREFIX.*"
echo "  Chromosome Dictionary: $CHROM_DICT"
echo "  Simplified GFF: $SIMPLIFIED_GFF"
echo "  Configuration: $CONFIG_FILE"
echo ""
echo "To use in your pipeline:"
echo "  1. Update GENOME_INDEX in process_scrna.sh:"
echo "     GENOME_INDEX=\"$PWD/$INDEX_PREFIX\""
echo ""
echo "  2. For gene mapping, use chromosome dictionary:"
echo "     $PWD/$CHROM_DICT"
echo ""
echo "Log files are available in the logs directory."
echo "==========================================================="

echo "[$(date)] Genome preparation pipeline completed successfully!"
