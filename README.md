This pipeline processes HIVE-style scRNA-seq data through a complete workflow including read merging, alignment, BAM processing, and gene quantification. It's optimized for high-throughput processing with multi-threading support and memory-efficient streaming algorithms.

Workflow:
Raw FASTQ Files (R1/R2)
        ↓
1. Merge R1/R2 + Extract Barcodes/UMIs
        ↓
2. Align to Reference Genome (HISAT2)
        ↓
3. BAM Processing & Sorting
        ↓
4. Stream Parse BAM Files
        ↓
5. Map Reads to Genes
        ↓
Gene × Cell Expression Matrix
