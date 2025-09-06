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
