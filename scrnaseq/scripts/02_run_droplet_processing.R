run_droplet_processing <- function(fpath, out_save = FALSE, plot_save = FALSE, force = FALSE) {
  message("--- Droplet processing ---")
  fname <- basename(tools::file_path_sans_ext(fpath))
  outfile <- paste0("scrnaseq/", fname, "_droplet.rds")
  
  if (file.exists(outfile)) {
    if (force) {
      break
    } else {
      message(paste0("File: '", outfile, "' already exists. Skipping."))
      return(NULL)
    }
  }
  message(paste0("Opening dataset '", fname, "'..."))
  sce <- DropletUtils::read10xCounts(paste0(fpath, "/raw_feature_bc_matrix"), col.names = TRUE, row.names = "symbol")
  bc_ranks <- DropletUtils::barcodeRanks(counts(sce))
  if (plot_save) {
    svg(filename = paste0("scrnaseq/plots/", fname, "_barcode_ranks.svg"))
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
    readr::write_rds(sce, outfile)
  }
  return(sce)
}