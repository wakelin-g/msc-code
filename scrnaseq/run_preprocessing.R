set.seed(42)

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  folders <- fs::dir_ls(args[[1]])
} else if (Sys.getenv("RSTUDIO") == "1") {
  folders <- fs::dir_ls("scrnaseq/outs/")
} else {
  setwd("/Users/griffen/Documents/thesis-code/")
  folders <- fs::dir_ls("scrnaseq/outs/")
}

suppressPackageStartupMessages(c(
  library(scater),
  library(scran),
  library(scuttle),
  library(SingleCellExperiment)
))

for (f in fs::dir_ls("scrnaseq/scripts")) {
  source(f)
}

bpp <- BiocParallel::MulticoreParam(8)

run_all <- function(sce, out_save = FALSE, plot_save = FALSE, force = FALSE) {
  run_soupx(sce, out_save, plot_save, force)
  run_droplet_processing(sce, out_save, plot_save, force)
  run_quickqc(sce, out_save, plot_save, force)
  run_normalization_dimred(sce, out_save, plot_save, force)
  run_clustering(sce, out_save, plot_save, force)
}

for (i in folders) {
  run_all(i, TRUE, TRUE, FALSE)
}
