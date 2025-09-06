create_S.O<- function(sample, group_name = ""){
  dataTran <- read.table(CMfile[grep(sample, CMfile)],
                         header=1,
                         row.names=1,
                         check.names=FALSE , 
                         sep="\t", 
                         quote = "", 
                         stringsAsFactors=0) 
  dataReads <- read.table(readsFile[grep(sample, readsFile)],
                          header=1,
                          row.names=1)
  table(rownames(dataReads) == colnames(dataTran))
  reads=dataReads[rownames(dataReads) %in% colnames(dataTran),]
  obj <- CreateSeuratObject(counts=dataTran, project = sample, min.cells = 10, min.features = 100)
  
  # Calculate metadata info
  logtotreads <- log10(as.matrix(reads)[,1])
  readAll <- dataReads[,1]
  readMap <- reads[,2]
  readMapPct <- readMap/readAll
  readUnmap <- readAll-readMap
  readUnmapPct <- readUnmap/readAll
  readExon <- reads[,3]
  readExonPct <- readExon/readAll
  readExonMapPct <- readExon/readMap
  names(readExon) <- colnames(obj)
  obj <- AddMetaData(obj, readExon, "ExonReads")
  names(readExonMapPct) <- colnames(obj)
  names(readExonPct) <- colnames(obj)
  names(readMap) <- colnames(obj)
  names(readAll) <- colnames(obj)
  
  # Add metadata to object
  obj <- AddMetaData(obj, logtotreads, "Log.TotReads") 
  obj <- AddMetaData(obj, readExonMapPct, "ExonvMapped")
  obj <- AddMetaData(obj, readExonPct, "ExonvTotal")
  obj <- AddMetaData(obj, readMap, "reads.mapped")
  obj <- AddMetaData(obj, readAll, "reads.Total")
  
  # Calculate complexity and add to metadata
  comp=obj@meta.data$ExonReads/obj@meta.data$nCount_RNA
  names(comp)<-rownames(obj@meta.data)
  obj <- AddMetaData(obj, comp, "Complexity")
  
  # Add group and groupSamp
  obj$group <- group_name
  obj$groupSamp <- paste(group_name, sample, sep = "_")
  
  return(obj)
}

create_S.O_merged<- function(sampleslist, groupslist){
  S.O.list <- purrr::map2(sampleslist, groupslist, 
                          ~ create_S.O(sample = .x, group_name = .y))
  if(length(S.O.list) > 1) {
    S.O.RNA <- merge(x = S.O.list[[1]], 
                     y = S.O.list[2:length(S.O.list)], 
                     add.cell.ids = sampleslist)
  } else {
    S.O.RNA <- S.O.list[[1]]
  }
  return(S.O.RNA)
}

clean_S.O <- function(S.O) {
  VlnPlot(S.O, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
  qfeature <- quantile(S.O$nFeature_RNA, probs = c(0.1, 0.98))
  qcount <- quantile(S.O$nCount_RNA, probs = c(0.1, 0.98))
  
  # Print thresholds
  cat("Thresholds used:\n")
  cat("nFeature_RNA:", qfeature[1], "to", qfeature[2], "\n")
  cat("nCount_RNA:", qcount[1], "to", qcount[2], "\n")

  # Filter cells
  S.O.subset <- subset(S.O, subset = nFeature_RNA > qfeature[1] & nFeature_RNA < qfeature[2] & 
                         nCount_RNA > qcount[1] & nCount_RNA < qcount[2])
  
  cat("\nCells before filtering:", ncol(S.O), "\n")
  cat("Cells after filtering:", ncol(S.O.subset), "\n")
  
  return(S.O.subset)
}

findNumberVariableGenes <- function(S.O, quantile_cutoff = 0.75) {
  S.O <- FindVariableFeatures(S.O, 
                              selection.method = "vst",
                              nfeatures = 6000)
  var_data <- HVFInfo(S.O)
  var_data <- var_data[order(var_data$variance.standardized, decreasing = TRUE), ]
  cutoff_value <- quantile(var_data$variance.standardized, quantile_cutoff)
  genes_above_cutoff <- sum(var_data$variance.standardized >= cutoff_value)
  return(genes_above_cutoff)
}

prep_S.O <- function(S.O, res = 0.5, find.ngenes = TRUE, cutoff = 0.75){
  set.seed(100)
  S.O <- NormalizeData(S.O, normalization.method = "LogNormalize", scale.factor = 10000)
  if (find.ngenes == TRUE){
    ngenes <- findNumberVariableGenes(S.O, cutoff)
  }
  else{
    ngenes = 5300
  }
  S.O <- FindVariableFeatures(S.O, selection.method = "vst", nfeatures = ngenes)
  all.genes <- rownames(S.O)
  S.O <- ScaleData(S.O, features = all.genes)
  S.O <- RunPCA(S.O, features = VariableFeatures(object = S.O), verbose = FALSE)
  S.O <- FindNeighbors(S.O, dims = 1:20, reduction = 'pca')
  S.O <- FindClusters(S.O, resolution = res)
  S.O <- RunUMAP(S.O, 
                 dims = 1:20,
                 n.components = 3,
                 n.neighbors = 25,     # Reduced to emphasize local structure
                 min_dist = 0.4,       # Increased to push points apart more
                 spread = 2.5,           # Increased to emphasize gaps
                 repulsion.strength = 2,  # Stronger repulsion to create hole
                 metric = 'correlation')  # Help push clusters apart
  return(S.O)
}

create_3DPlot <- function(S.O){
  embeddings <- Embeddings(S.O, "umap")
  plot_data <- data.frame(
    UMAP_1 = embeddings[,1],
    UMAP_2 = embeddings[,2],
    UMAP_3 = embeddings[,3],
    cluster = Idents(S.O)
  )
  p <- plot_ly(plot_data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~cluster, 
          type = "scatter3d", mode = "markers", marker = list(size = 2))
  p <- p %>% layout(
    legend = list(itemsizing = "constant", itemwidth = 30)
  )
  p
}

plotUMAP <- function(S.O, x, y) {
  umap_data <- Embeddings(S.O, "umap")
  plot_data <- data.frame(
    UMAP_1 = umap_data[,x],  
    UMAP_2 = umap_data[,y],  
    cluster = as.factor(Idents(S.O))  
  )
  
  ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
    geom_point(alpha = 0.7) +
    theme_classic() +
    labs(xaxis = list(title = paste("UMAP Dimension ", x)),
         yaxis = list(title = paste("UMAP Dimension ", y)))
}

topFeaturePlot <- function(S.O, gene, top_percent){
  library(ggplot2)
  # Extract UMAP coordinates and gene expression
  umap_coords <- Embeddings(S.O, reduction = "umap")
  gene_id <- gene  # Adjust to match your exact gene ID
  gene_expr <- FetchData(S.O, vars = gene_id)[, 1]
  
  # Create a data frame with coordinates and expression
  plot_data <- data.frame(
    UMAP_2 = umap_coords[, 2],
    UMAP_3 = umap_coords[, 3],
    Expression = gene_expr
  )
  
  # Keep only cells with expression above a threshold (e.g., top 10%)
  percent <- 1 - top_percent
  threshold <- quantile(gene_expr[gene_expr > 0], percent)  # Only top 10% of non-zero expressing cells
  filtered_data <- plot_data[plot_data$Expression >= threshold, ]
  
  # Show all cells but highlight high expressors
  ggplot() +
    # First layer: all cells in gray
    geom_point(data = plot_data, aes(x = UMAP_2, y = UMAP_3), 
               color = "lightgray", size = 0.1, alpha = 0.3) +
    # Second layer: only high-expressing cells colored by expression
    geom_point(data = filtered_data, aes(x = UMAP_2, y = UMAP_3, color = Expression), 
               size = 1) +
    scale_color_viridis_c() +
    theme_classic() +
    labs(title = paste0(gene_id, " (highlighted: top ", round(nrow(filtered_data)/nrow(plot_data)*100), "% expressing cells)"))
}

add_gene_descriptions_from_gff_fast <- function(all_markers, gff_path, gene_column = "gene") {
  # Read GFF manually (ignoring comment lines)
  gff <- read.delim(gff_path, header = FALSE, comment.char = "#", sep = "\t", stringsAsFactors = FALSE)
  colnames(gff)[9] <- "attributes"
  
  # Filter only protein_coding_gene
  gff_genes <- gff[gff$V3 == "protein_coding_gene", ]
  
  # Extract ID and description from attributes field
  extract_field <- function(attr, field) {
    pattern <- paste0(field, "=([^;]*)")
    match <- regmatches(attr, regexpr(pattern, attr))
    sub(paste0(field, "="), "", match)
  }
  
  ids <- extract_field(gff_genes$attributes, "ID")
  descriptions <- extract_field(gff_genes$attributes, "description")
  
  # Convert ID from underscore to dash to match all_markers$gene
  ids <- gsub("_", "-", ids)
  
  # Create lookup table
  id_to_description <- setNames(descriptions, ids)
  
  # Add description to all_markers
  all_markers$gene_name <- id_to_description[as.character(all_markers[[gene_column]])]
  
  return(all_markers)
}

### creates and returns pseudotime expression matrix. To save pseudotime in object, path should be passed as overwrite
### can only be used for simple pseudotime curves (no bifurcation)
pseudotimeExpressionSimple<-function(seurat_object, start_cluster='ring', overwrite=NULL){
  data<-seurat_object
  
  DefaultAssay(data)<-'RNA'
  sling_data <- slingshot(Embeddings(data, "pca"), clusterLabels = data$stage, start.clus = start_cluster)
  sling_data <- SlingshotDataSet(sling_data)
  plot3d.SlingshotDataSet(sling_data)
  cell_id <- as.data.frame(slingBranchID(sling_data))
  
  pseudotime2 <- as.data.frame(slingPseudotime(sling_data))
  pseudotime2$cell <- rownames(pseudotime)
  pseudotime2$stage <- NA
  for (i in 1:ncol(sling_data@clusterLabels)){
    pseudotime[which(pseudotime$cell %in% rownames(sling_data@clusterLabels)[which(sling_data@clusterLabels[, i] == 1)]), 
               'stage'] <- colnames(sling_data@clusterLabels)[i]
  }
  
  if (is.null(overwrite) == FALSE){
    ind <- colnames(data) %in% pseudotime$cell
    idc <- pseudotime$Lineage1
    data@meta.data$Pseudotime <- ifelse(colnames(data) %in% pseudotime$cell, idc, NA)
    saveRDS(data, overwrite)
  }
  
  gene_exp <- as.data.frame(data@assays$RNA@counts)
  ind <- which(pseudotime$cell %in% rownames(cell_id))
  idc_pseudotime <- pseudotime[ind, c("Lineage1", "cell", 'stage')]
  temp <- idc_pseudotime[order(idc_pseudotime$Lineage1), ]
  ordered_cells <- rownames(temp)
  idc_pseudotime <- as.data.frame(t(idc_pseudotime))
  ind2 <- which(colnames(gene_exp) %in% rownames(cell_id))
  idc_pseudotime <- rbind(idc_pseudotime, gene_exp[, ind2])
  idc_pseudotime <- idc_pseudotime[, ordered_cells]
  
  return(idc_pseudotime)
}



plotPseudotimeHeatmap<-function(pseudotime_expression, genes_interest, list_colors, list_stages=c('ring','trophozoite','schizont')){
  cells<-unlist(as.vector(pseudotime_expression['stage',]))
  pseudo_exp<-pseudotime_expression[c(4:nrow(pseudotime_expression)),]
  genes<-rownames(pseudo_exp)
  pseudo_exp <- as.data.frame(sapply(pseudo_exp, as.numeric))
  rownames(pseudo_exp)<-genes
  pseudo_exp<-pseudo_exp[]+1
  pseudo_exp_log<-log(pseudo_exp)
  ind<-which(rownames(pseudo_exp) %in% genes_interest)
  my_cols<-list_colors[match(cells, list_stages)]
  heat.map<-heatmap(data.matrix(pseudo_exp_log[ind,]), Rowv = NA, Colv = NA, scale="column", ColSideColors =my_cols)
  
  return(heat.map)
}


### creates a matrix of predicted expression values on pseudotime using splines
predictedExpressionMatrix<-function(seurat_object, pseudotime_expression, genes_interest){
  data<-seurat_object
  pseudo_exp_file<-pseudotime_expression
  goi<-genes_interest  
  
  exp_mat<-as.matrix(data@assays$RNA$data)
  exp_mat_filter<-exp_mat[goi,]
  cells<-colnames(pseudo_exp_file)[which(colnames(pseudo_exp_file) %in% colnames(exp_mat_filter))]
  exp_mat_filter<-exp_mat_filter[,cells]
  
  max_time<-as.numeric(pseudo_exp_file[1,ncol(pseudo_exp_file)])
  pseudo_exp<-pseudo_exp_file[c(4:nrow(pseudo_exp_file)),]
  genes<-rownames(pseudo_exp)
  pseudo_exp <- as.data.frame(sapply(pseudo_exp, as.numeric))
  rownames(pseudo_exp)<-genes
  
  pseudotime<-as.numeric(pseudo_exp_file[1,which(colnames(pseudo_exp_file) %in% colnames(exp_mat_filter))])
  scaled_pseudo<-pseudotime/max_time
  
  first_gene<-exp_mat_filter[goi[1],]
  first_gene_spline<-data.frame(scaled_pseudo, first_gene)
  spline<-smooth.spline(first_gene_spline$scaled_pseudo, first_gene_spline$first_gene)
  prediction <- predict(spline, seq(0, 1, by = 0.01))
  exp_mat_filter_pred<-as.data.frame(prediction$y)
  
  for (i in 2:nrow(exp_mat_filter)){
    new_gene<-exp_mat_filter[goi[i],]
    new_gene_spline<-data.frame(scaled_pseudo, new_gene)
    spline<-smooth.spline(new_gene_spline$scaled_pseudo, new_gene_spline$new_gene)
    prediction <- predict(spline, seq(0, 1, by = 0.01))
    exp_mat_filter_pred<-cbind(exp_mat_filter_pred, prediction$y)
  }
  
  colnames(exp_mat_filter_pred)<-goi
  exp_mat_filter_pred<-t(exp_mat_filter_pred)
  colnames(exp_mat_filter_pred)<-seq(0,1,0.01)
  
  return(exp_mat_filter_pred)
}


### plots predicted expression on pseudotime with colored bar for stages
plotExpressionPseudotime<-function(seurat_object, predicted_expression, colors_list, 
                                   stages_list=c('schizont','trophozoite','ring'), plot_title=NULL){
  
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(cowplot)
  
  data<-seurat_object
  exp_mat_filter_pred<-pseudotime_expression
  
  stages_table<-table(data$stage)
  color_bar<-data.frame(x=1,y=as.numeric(stages_table), 
                        z=names(stages_table))
  color_bar$z<-factor(color_bar$z, levels=stages_list)
  
  colored_plot<-ggplot(color_bar, aes(x=x, y=y, fill=z)) + 
    geom_bar(position="stack", stat="identity", show.legend = FALSE) +
    scale_fill_manual(values=colors_list) +
    coord_flip()+
    theme_void()
  
  final_exp<-as.data.frame(t(exp_mat_filter_pred))
  final_exp_long<-melt(t(as.data.frame(exp_mat_filter_pred)))
  colnames(final_exp_long)<-c('Pseudotime', 'Gene', 'Expression')
  data_ends<-final_exp_long %>% filter(Pseudotime==1.00)
  
  exp_plot<-ggplot(final_exp_long, aes(x=Pseudotime, y=Expression, group=Gene)) + 
    geom_line(aes(color=Gene),linewidth=2) +
    theme_minimal() +
    ggtitle(plot_title)
  
  final_plot<-plot_grid(exp_plot, colored_plot, ncol = 1, rel_heights = c(10, .5))
  
  return(final_plot)
  
}
