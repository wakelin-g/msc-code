run_quickqc <- function(fpath, out_save = FALSE, plot_save = FALSE, force = FALSE) {
  message("--- Quality control ---")
  fname <- basename(tools::file_path_sans_ext(fpath))
  infile <- paste0("scrnaseq/", fname, "_droplet.rds")
  outfile <- paste0("scrnaseq/", fname, "_droplet_qc.rds")
  
  if (file.exists(outfile)) {
    if (force) {
      break
    } else {
      message("File: '", outfile, "' already exists. Skipping.")
      return(NULL)
    }
  }
  
  sce <- readr::read_rds(infile)
  qc_stats <- scuttle::perCellQCMetrics(sce, subsets = list(Mt = grepl("mt-", rownames(sce))))
  qc_filters <- scuttle::perCellQCFilters(qc_stats, nmads = 3)
  colData(sce) <- cbind(colData(sce), qc_stats)
  colData(sce) <- cbind(colData(sce), qc_filters)
  colData(sce)$cell_doublet_scores <- scDblFinder::computeDoubletDensity(sce)
  
  if (plot_save) {
    if (!dir.exists(paste0("scrnaseq/plots/", fname, "/"))) {
      dir.create(paste0("scrnaseq/plots/", fname))
    }
    for (qc_stat in c("subsets_Mt_percent", "sum", "detected")) {
      p <- scater::plotColData(sce, y = qc_stat)
      ggplot2::ggsave(plot = p, filename = paste0("scrnaseq/plots/", fname, "/", qc_stat, ".svg"))
    }
  }
  
  readr::write_rds(sce, outfile)
  return(sce)
}