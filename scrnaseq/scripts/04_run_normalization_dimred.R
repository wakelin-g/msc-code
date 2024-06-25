run_normalization_dimred <- function(fpath, out_save = FALSE, plot_save = FALSE, force = FALSE) {
  message("--- Dimensionality reduction ---")
  fname <- basename(tools::file_path_sans_ext(fpath))
  infile <- paste0("scrnaseq/", fname, "_droplet_qc.rds")
  outfile <- paste0("scrnaseq/", fname, "_droplet_qc_norm_dimred.rds")
 
  if (file.exists(outfile))  {
    if (force) {
      break
    } else {
      message(paste0("File: '", outfile, "' already exists. Skipping."))
      return(NULL)
    }
  }
  
  sce <- readr::read_rds(infile)
  
  clusters <- scran::quickCluster(sce)
  sce <- scran::computeSumFactors(sce, clusters = clusters)
  sce <- scuttle::logNormCounts(sce)
  
  gene_var <- scran::modelGeneVarByPoisson(sce)
  gene_hvgs <- scran::getTopHVGs(sce, prop = 0.1)
  
  sce <- scran::denoisePCA(sce, subset.row = gene_hvgs, technical = gene_var)
  sce <- scater::runTSNE(sce, dimred = "PCA")
  sce <- scater::runUMAP(sce, dimred = "PCA")
  
  if (plot_save) {
    svg(filename = paste0("scrnaseq/plots/", fname, "/gene_var.svg"))
    plot(gene_var$mean, gene_var$total, pch = 16, cex = 0.5, xlab = "Mean of log-expression", ylab = "Variance of log-expression")
    curfit <- metadata(gene_var)
    curve(curfit$trend(x), col = "red", add = TRUE, lwd = 2)
    dev.off()
    rm(curfit)
  }
  
  if (plot_save) {
    if (!dir.exists(paste0("scrnaseq/plots/", fname, "/"))) {
      dir.create(paste0("scrnaseq/plots/", fname))
    }
    for (qc_stat in c("subsets_Mt_percent", "sum", "detected")) {
      p <- scater::plotColData(sce, y = qc_stat)
      ggplot2::ggsave(plot = p, filename = paste0("scrnaseq/plots/", fname, "/", qc_stat, ".svg"))
    }
  }
  
  if (out_save) {
    readr::write_rds(sce, outfile)
  }
  return(sce)
}