#' runEmbedding
#'
#' @param req type = (tsne,umap,pca)
#'            umap config = (minimumDistance, distanceMetric)
#'            tsne config = (perplexity,learningRate)
#' @param data Object of class Seurat
#'
#' @return x,y positions for each cell, sorted by cellId with NA values in filtered cells.
#' @export
#'
#' @example As
runEmbedding <- function(req, data) {
  type <- req$body$type
  config <- req$body$config
  pca_nPCs <- 30

  # To run embedding, we need to set the reduction.
  if ("active.reduction" %in% names(data@misc)) {
    active.reduction <- data@misc[["active.reduction"]]
  } else {
    active.reduction <- "pca"
  }

  # The slot numPCs is set in dataIntegration with the selectd PCA by the user.
  if ("numPCs" %in% names(data@misc)) {
    pca_nPCs <- data@misc[["numPCs"]]
  }

  message("Active reduction --> ", active.reduction)
  message("Active numPCs --> ", pca_nPCs)
  message("Number of cells/sample:")
  table(data$samples)

  if (type == "pca") {
    # Leaving this here to add parameters in the future. Won't leave uncommented to avoid recalculating PCA
    # RunPCA(data, npcs = 50, features = VariableFeatures(object=data), verbose=FALSE)
    df_embedding <- Seurat::Embeddings(data, reduction = type)[, 1:2]
  } else if (type == "tsne") {
    data <- Seurat::RunTSNE(data,
      reduction = active.reduction,
      seed.use = 1,
      dims = 1:pca_nPCs,
      perplexity = config$perplexity,
      learning.rate = config$learningRate
    )
    df_embedding <- Seurat::Embeddings(data, reduction = type)
  } else if (type == "umap") {
    data <- Seurat::RunUMAP(data,
      seed.use = 42,
      reduction = active.reduction,
      dims = 1:pca_nPCs,
      verbose = FALSE,
      min.dist = config$minimumDistance,
      metric = config$distanceMetric,
      umap.method = "umap-learn"
    )

    df_embedding <- Seurat::Embeddings(data, reduction = type)
  }

  # Order embedding by cells id in ascending form
  df_embedding <- as.data.frame(df_embedding)
  df_embedding$cells_id <- data@meta.data$cells_id
  df_embedding <- df_embedding[order(df_embedding$cells_id), ]
  df_embedding <- df_embedding %>%
    tidyr::complete(cells_id = seq(0, max(data@meta.data$cells_id))) %>%
    select(-cells_id)

  map2_fun <- function(x, y) {
    if (is.na(x)) {
      return(NULL)
    } else {
      return(c(x, y))
    }
  }
  res <- purrr::map2(df_embedding[[1]], df_embedding[[2]], map2_fun)
  return(res)
}
