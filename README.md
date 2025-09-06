This pipeline processes HIVE-style scRNA-seq data through a complete workflow including read merging, alignment, BAM processing, and gene quantification (Seurat objet creation is optional). It's optimized for high-throughput processing with multi-threading support and memory-efficient streaming algorithms. The run is done in four main steps, each handling specific aspects of the analysis workflow.

General workflow:
Raw FASTQ Files (R1/R2) → Merge R1/R2 + Extract Barcodes/UMIs → Align to Reference Genome (HISAT2) → BAM Processing & Sorting → Stream Parse BAM Files →  Map Reads to Genes → 
Gene × Cell Expression Matrix → Create Seurat object (optional)

Steps to run the pipeline:
### 1. Prepare Reference Genome (once per organism)
	# Creates a folder named GenomeName, should be run inside Genomes folder
		./prepare_genome.sh \
  			--fasta https://example.com/genome.fa.gz \
	  		--gff https://example.com/annotations.gff3.gz \
  			--output GenomeName 

### 2. Process Each Lane: 
	# Run for each lane separately
		./process_scrna.sh -1 lane1_R1.fastq.gz -2 lane1_R2.fastq.gz > process.log 2>&1 &
		./process_scrna.sh -1 lane2_R1.fastq.gz -2 lane2_R2.fastq.gz > process.log 2>&1 &
		# ... repeat for all lanes

### 3. Merge Multi-Lane Results
	# Merges sorted BAM file (outputed by process_scrna) and parses it creating a merged_duplicate_reads.txt table
		./merge_lanes_and_parse.sh lane1/ lane2/ lane3/ > process.log 2>&1 &

### 4. Create Expression Matrix
	# Gets merged_duplicate_reads.txt and creates a Gene x Cell matrix to be exported to R
		 python -u ./create_count_matrix.py duplicate_reads.txt chromosome_dict.pkl -o expression_matrix.txt -m 5000 > mapping_output.log 2>&1 &
