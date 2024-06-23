folders <- fs::dir_ls("outs/")
soupx_outs_dir <- "outs_soupx/"

for (i in folders) {
  fname <- strsplit(i, "/")[[1]][[2]]
  sc <- SoupX::load10X(i)
  sc <- SoupX::autoEstCont(sc)
  out <- SoupX::adjustCounts(sc)
  Matrix::writeMM(out, paste0(soupx_outs_dir, fname, "_soupx_corrected.mtx"))
}
