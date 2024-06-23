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
  fname = basename(fpath)
  message(paste0("Opening dataset '", fname, "'..."))
  sce <- DropletUtils::read10xCounts(paste0(fpath, "/raw_feature_bc_matrix"), col.names = TRUE, row.names = "symbol")
  bc_ranks <- DropletUtils::barcodeRanks(counts(sce))
  if (plot_save) {
    svg(filename = paste0("outs_droplets/", fname, "_barcode_ranks.svg"))
    bc_ranks_unique <- !duplicated(bc_ranks$rank)
    plot(bc_ranks$rank[bc_ranks_unique], bc_ranks$total[bc_ranks_unique], log = "xy", xlab = "Rank", ylab = "Total UMI Count", cex.lab = 1.2, main = fname)
    abline(h = metadata(bc_ranks)$inflection, col = "darkgreen", lty = 2)
    abline(h = metadata(bc_ranks)$knee, col = "skyblue", lty = 2)
    legend("bottomleft", legend = c("Inflection", "Knee"), col = c("darkgreen", "skyblue"), lty = 2, cex = 1.2)
    dev.off()
  }
  
  empty_out <- DropletUtils::emptyDrops(counts(sce), lower = 50, test.ambient = TRUE, BPPARAM = bpp)
  sce <- sce[,which(empty_out$FDR <= 0.001)]
  
  if (out_save) {
    readr::write_rds(sce, paste0(fname, "_droplet.rds"))
  }
  return(sce)
}

run_quickqc <- function(fpath_sce) {
  sce <- readr::read_rds(fpath_sce)
  stats <- scater::perCellQCMetrics(sce, subsets = list(Mt = grepl("mt-", rownames(sce))))
  
  mt_high <- scater::isOutlier(stats$subsets_Mt_percent, type = "higher")
  
  sce_out <- sce[,!discard]
  
  readr::write_rds(sce_out, paste0(strsplit(fpath_sce, ".")[[1]][[1]], "_qc.rds"))
}

run_quickclustering <- function(fpath_sce) {
  sce <- readr::read_rds(fpath_sce)
  
  gene_hvgs <- scran::getTopHVGs(sce, prop = 0.1)
  gene_var <- scran::modelGeneVarByPoisson(sce)
  
  sce <- scran::denoisePCA(sce, subset.row=gene_hvgs, technical=gene_var)
  sce <- scater::runUMAP(sce, use.dimred = "PCA")
  
  
  readr::write_rds(sce_out, paste0(strsplit(fpath_sce, ".")[[1]][[1]], "_clustered.rds"))
}

for (i in folders) {
  fname <- strsplit(i, "/")[[1]][[2]]

  run_soupx(i, out_save = TRUE, plot_save = TRUE)
  
  #sce <- DropletUtils::read10xCounts(paste0(i, "/filtered_feature_bc_matrix"))
  #rownames(sce) <- SummarizedExperiment::rowData(sce)$Symbol
  #SummarizedExperiment::assay(sce, "raw") <- SingleCellExperiment::counts(sce)
  #SummarizedExperiment::assay(sce, "counts", withDimnames = FALSE) <- soup

  #qc_stats <- scuttle::perCellQCMetrics(sce, subsets = list(Mt = grepl("mt-", rownames(sce))))
  #qc_high <- scuttle::isOutlier(qc_stats$subsets_Mt_percent, type = "higher", nmads = 3)
  #discard <- scuttle::quickPerCellQC(qc_stats, percent_subsets = c("subsets_Mt_percent"))
  #SummarizedExperiment::colData(sce) <- cbind(SummarizedExperiment::colData(sce), qc_stats)
  # sce <- sce[, !discard$discard]
  #
  #clusters <- scran::quickCluster(sce)
  #sce <- scran::computeSumFactors(sce, cluster = clusters)
  #sce <- scuttle::logNormCounts(sce)

  #gene_var <- scran::modelGeneVarByPoisson(sce, BPPARAM = bpp)
  #hvgs_top <- scran::getTopHVGs(sce, prop = 0.1)

  #sce <- scran::denoisePCA(sce, technical = gene_var, subset.row = hvgs_top)
  #sce <- scater::runTSNE(sce, dimred = "PCA")
  #sce <- scater::runUMAP(sce, dimred = "PCA")
  #sce <- scater::runMDS(sce, dimred = "PCA")

  #hashed.doublets <- scran::recoverDoublets(sce, use.dimred = "PCA", doublets = sce$Doublet)

  #snn <- scran::buildSNNGraph(sce, use.dimred = "PCA", k = 25)
  #SingleCellExperiment::colLabels(sce) <- factor(igraph::cluster_walktrap(snn)$membership)

  #readr::write_rds(sce, paste0(fname, "_sce.rds"))
}


















