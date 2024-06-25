run_clustering <- function(fpath, out_save = FALSE, plot_save = FALSE, force = FALSE) {
  message("--- Clustering ---")
  fname <- basename(tools::file_path_sans_ext(fpath))
  infile <- paste0("scrnaseq/", fname, "_droplet_qc_norm_dimred.rds")
  outfile <- paste0("scrnaseq/", fname, "_droplet_qc_norm_dimred_clustered.rds")
 
  if (file.exists(outfile))  {
    if (force) {
      break
    } else {
      message(paste0("File: '", outfile, "' already exists. Skipping."))
      return(NULL)
    }
  }
  
  sce <- readr::read_rds(infile)
  
  snn <- scran::buildSNNGraph(sce, use.dimred = "PCA", k = 25)
  colData(sce)$walktrap <- factor(igraph::cluster_walktrap(snn, steps = 4)$membership)
  colData(sce)$leiden <- factor(igraph::cluster_leiden(snn, resolution_parameter = 0.5)$membership)
  
  if (plot_save) {
    for (cluster_type in c("leiden", "walktrap")) {
      for (plot_type in c("PCA", "TSNE", "UMAP")) {
        do.call(paste0("plot", plot_type), list(object = sce, color_by = cluster_type))
        ggplot2::ggsave(filename = paste0("scrnaseq/plots/", fname, "/", plot_type, "_", cluster_type, ".svg"))
      }
    }
  }
  
  #markers.leiden <- scran::findMarkers(sce, groups = sce$leiden, pval.type = "all", direction = "up")
  #markers.walktrap <- scran::findMarkers(sce, groups = sce$walktrap, pval.type = "all", direction = "up")
  markers.idf.leiden <- SoupX::quickMarkers(counts(sce), clusters = sce$leiden, N = 10, FDR = 0.01)
  markers.idf.walktrap <- SoupX::quickMarkers(counts(sce), clusters = sce$walktrap, N = 10, FDR = 0.01)
  
  if (out_save) {
    readr::write_rds(sce, outfile)
    readr::write_csv(markers.idf.leiden, paste0(tools::file_path_sans_ext(outfile), "_markers_leiden.csv"))
    readr::write_csv(markers.idf.walktrap, paste0(tools::file_path_sans_ext(outfile), "_markers_walktrap.csv"))
  }
  
  return(sce)
}