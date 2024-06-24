set.seed(42)

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  folders <- fs::dir_ls(args[[1]])
} else {
  folders <- fs::dir_ls("outs/")
}

suppressPackageStartupMessages(library(SingleCellExperiment))

bpp <- BiocParallel::MulticoreParam(8)

run_soupx <- function(fpath, out_save = TRUE, plot_save = FALSE) {
  message("--- Ambient RNA estimation ---")
  fname = basename(fpath)
  message(paste0("Opening dataset '", fname, "'..."))
  soup <- SoupX::load10X(fpath, verbose = FALSE)
  if (plot_save) {
    if (!dir.exists("outs_soupx")) {
      dir.create("outs_soupx")
    }
    message("Estimating contamination; saving SoupX plot.")
    svg(filename = paste0("outs_soupx/", fname, "_soupx_est_contamination.svg"))
    soup <- SoupX::autoEstCont(soup, doPlot = TRUE, verbose = FALSE)
    dev.off()
  } else {
    message("Estimating contamination; not saving SoupX plot.")
    soup <- SoupX::autoEstCont(soup, doPlot = FALSE, verbose = FALSE)
  }
  message("Adjusting counts based on estimated contamination.")
  soup_out <- SoupX::adjustCounts(soup, verbose = 0)
  if (out_save) {
    message(paste0("Saving SoupX outs to 'outs_soupx/", fname, "_soupx.mtx'..."))
    if (!dir.exists("outs_soupx")) {
      dir.create("outs_soupx")
    }
    Matrix::writeMM(soup_out, paste0("outs_soupx/", fname, "_soupx.mtx"))
  }
  return(soup_out)
}

run_droplet_processing <- function(fpath, out_save = TRUE, plot_save = FALSE) {
  message("--- Droplet processing ---")
  fname = basename(fpath)
  message(paste0("Opening dataset '", fname, "'..."))
  sce <- DropletUtils::read10xCounts(paste0(fpath, "/raw_feature_bc_matrix"), col.names = TRUE, row.names = "symbol")
  bc_ranks <- DropletUtils::barcodeRanks(counts(sce))
  if (plot_save) {
    if (!dir.exists("outs_droplets")) {
      dir.create("outs_droplets")
    }
    svg(filename = paste0("plots/", fname, "_barcode_ranks.svg"))
    bc_ranks_unique <- !duplicated(bc_ranks$rank)
    plot(bc_ranks$rank[bc_ranks_unique], bc_ranks$total[bc_ranks_unique], log = "xy", xlab = "Rank", ylab = "Total UMI Count", cex.lab = 1.2, main = fname)
    abline(h = metadata(bc_ranks)$inflection, col = "darkgreen", lty = 2)
    abline(h = metadata(bc_ranks)$knee, col = "skyblue", lty = 2)
    legend("bottomleft", legend = c("Inflection", "Knee"), col = c("darkgreen", "skyblue"), lty = 2, cex = 1.2)
    dev.off()
  }
  message("Calculating empty droplets...")
  empty_out <- DropletUtils::emptyDrops(counts(sce), lower = 50, test.ambient = TRUE, BPPARAM = bpp)
  sce <- sce[,which(empty_out$FDR <= 0.001)]
  
  if (out_save) {
    readr::write_rds(sce, paste0(fname, "_droplet.rds"))
  }
  return(sce)
}

run_quickqc <- function(fpath_sce) {
  message("--- Quality control ---")
  sce <- readr::read_rds(fpath_sce)
  stats <- scater::perCellQCMetrics(sce, subsets = list(Mt = grepl("mt-", rownames(sce))))
  
  mt_high <- scater::isOutlier(stats$subsets_Mt_percent, type = "higher")
  dbl_clusters <- scDblFinder::findDoubletClusters(sce)
  dbl_simulated <- scDblFinder::computeDoubletDensity(sce, subset.row=)
  
  
  
  sce_out <- sce[,!discard]
  
  readr::write_rds(sce_out, paste0(strsplit(fpath_sce, ".")[[1]][[1]], "_qc.rds"))
  return(sce_out)
}

run_normalization_dimred <- function(fpath_sce, out_save = TRUE, plot_save = FALSE) {
  message("--- Dimensionality reduction ---")
  sce <- readr::read_rds(fpath_sce)
  
  clusters <- scran::quickCluster(sce)
  sce <- scran::computeSumFactors(sce, clusters = clusters)
  sce <- scuttle::logNormCounts(sce)
  
  gene_var <- scran::modelGeneVarByPoisson(sce)
  gene_hvgs <- scran::getTopHVGs(sce, prop = 0.1)
  
  if (plot_save) {
    svg(filename = paste0("plots/", strsplit(fpath_sce, ".")[[1]][[1]], "_gene_var.svg"))
    plot(gene_var$mean, gene_var$total, pch = 16, cex = 0.5, xlab = "Mean of log-expression", ylab = "Variance of log-expression")
    curfit <- metadata(gene_var)
    curve(curfit$trend(x), col = "red", add = TRUE, lwd = 2)
    dev.off()
    rm(curfit)
  }
  
  sce <- scran::denoisePCA(sce, subset.row = gene_hvgs, technical = gene_var)
  sce <- scater::runTSNE(sce, dimred = "PCA")
  sce <- scater::runUMAP(sce, dimred = "PCA")
  
  if (out_save) {
    readr::write_rds(sce, paste0(strsplit(basename(fpath_sce), ".", fixed = TRUE)[[1]][[1]], "_norm_dimred.rds"))
  }
  return(sce)
}

run_clustering <- function(fpath_sce, out_save = TRUE, plot_save = FALSE) {
  message("--- Clustering ---")
  sce <- readr::read_rds(fpath_sce)
  
  snn <- scran::buildSNNGraph(sce, use.dimred = "PCA", k = 25)
  colData(sce)$walktrap <- factor(igraph::cluster_walktrap(snn, steps = 4)$membership)
  colData(sce)$leiden <- factor(igraph::cluster_leiden(snn, resolution_parameter = 0.5)$membership)
  
  if (plot_save) {
    svg(filename = paste0("plots/", strsplit(basename(fpath_sce), ".", fixed = TRUE)[[1]][[1]], "_pcaplot_walktrap.svg"))
    scater::plotPCA(sce, color_by="walktrap")
    dev.off()
    svg(filename = paste0("plots/", strsplit(basename(fpath_sce), ".", fixed = TRUE)[[1]][[1]], "_pcaplot_leiden.svg"))
    scater::plotPCA(sce, color_by="leiden")
    dev.off()
    svg(filename = paste0("plots/", strsplit(basename(fpath_sce), ".", fixed = TRUE)[[1]][[1]], "_tsneplot_walktrap.svg"))
    scater::plotTSNE(sce, color_by="walktrap")
    dev.off()
    svg(filename = paste0("plots/", strsplit(basename(fpath_sce), ".", fixed = TRUE)[[1]][[1]], "_tsneplot_leiden.svg"))
    scater::plotTSNE(sce, color_by="leiden")
    dev.off()
    svg(filename = paste0("plots/", strsplit(basename(fpath_sce), ".", fixed = TRUE)[[1]][[1]], "_umapplot_walktrap.svg"))
    scater::plotUMAP(sce, color_by="walktrap")
    dev.off()
    svg(filename = paste0("plots/", strsplit(basename(fpath_sce), ".", fixed = TRUE)[[1]][[1]], "_umapplot_leiden.svg"))
    scater::plotUMAP(sce, color_by="leiden")
    dev.off()
  }
  
  markers.leiden <- scran::findMarkers(sce, groups = sce$leiden, pval.type = "all", direction = "up")
  markers.walktrap <- scran::findMarkers(sce, groups = sce$walktrap, pval.type = "all", direction = "up")
  
  if (out_save) {
    readr::write_rds(sce, paste0(strsplit(basename(fpath_sce), ".", fixed = TRUE)[[1]][[1]], "_clustered.rds"))
    readr::write_rds(c(markers.walktrap, markers.leiden), paste0(strsplit(basename(fpath_sce), ".", fixed = TRUE)[[1]][[1]], "_clustered.rds"))
  }
  
  return(sce)
}

for (i in folders) {
  fname <- strsplit(i, "/")[[1]][[2]]

  run_soupx(i, out_save = TRUE, plot_save = TRUE)
  run_droplet_processing(i, out_save = TRUE, plot_save = TRUE)
  run_normalization_dimred(paste0(fname, "_droplet.rds"), plot_save = TRUE)
  run_clustering(paste0(fname, "_droplet_norm_dimred.rds"), plot_save = TRUE)
}
