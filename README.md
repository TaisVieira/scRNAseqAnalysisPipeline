This pipeline processes HIVE-style scRNA-seq data through a complete workflow including read merging, alignment, BAM processing, and gene quantification. It's optimized for high-throughput processing with multi-threading support and memory-efficient streaming algorithms.

Workflow:
Raw FASTQ Files (R1/R2) → Merge R1/R2 + Extract Barcodes/UMIs → Align to Reference Genome (HISAT2) → BAM Processing & Sorting → Stream Parse BAM Files →  Map Reads to Genes → 
Gene × Cell Expression Matrix
