tf_idf <- function(sce, clusters, expression_cutoff = 0.9) {
  x <- as(counts(sce), "dgCMatrix")
  w <- which(x@x > expression_cutoff)
  cluster_counts <- table(clusters)
  n_obs <- split(factor(rownames(x))[x@i[w]+1])
  
  
}