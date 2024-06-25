run_soupx <- function(fpath, out_save = TRUE, plot_save = FALSE, force = FALSE) {
  message("--- Ambient RNA estimation ---")
  fname <- basename(tools::file_path_sans_ext(fpath))
  outfile <- paste0("scrnaseq/outs_soupx/", fname, "_soupx.mtx")
  
  if (file.exists(outfile)) {
    if (force) {
      break
    } else {
      message("File: '", outfile, "' already exists. Skipping.")
      return(NULL)
    }
  }
  
  message(paste0("Opening dataset '", fname, "'..."))
  soup <- SoupX::load10X(fpath, verbose = FALSE)
  if (plot_save) {
    if (!dir.exists("scrnaseq/outs_soupx")) {
      dir.create("scrnaseq/outs_soupx")
    }
    message("Estimating contamination; saving SoupX plot.")
    svg(filename = paste0("scrnaseq/outs_soupx/", fname, "_soupx_est_contamination.svg"))
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
    if (!dir.exists("scrnaseq/outs_soupx")) {
      dir.create("scrnaseq/outs_soupx")
    }
    Matrix::writeMM(soup_out, file = outfile)
  }
  return(soup_out)
}