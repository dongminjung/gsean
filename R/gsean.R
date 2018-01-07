# GSEA with networks
gsean <- function(geneset, x, exprs, pseudo = 1, threshold = 0.99, nperm = 1000,
                  centrality = function(x) rowSums(abs(x)), weightParam = 1,
                  minSize = 1, maxSize = Inf, gseaParam = 1, nproc = 0,
                  BPPARAM = NULL, corParam = list(), tmax = 10, ...)
{
  if(!is.character(x) & !is.numeric(x))
    stop("list or statistic is not appropriate for ORA or GSEA")
  
  message("construct adjacency matrix ... ")
  adjacency <- do.call("exprs2adj", c(list(exprs, pseudo = pseudo), corParam))
  if(is.numeric(x))
  {
    message("GSEA ... ")
    centrality.gsea(geneset, x, adjacency, pseudo = pseudo, centrality = centrality,
                    weightParam = weightParam, nperm = nperm, minSize = minSize,
                    maxSize = maxSize, gseaParam = gseaParam, nproc = nproc,
                    BPPARAM = BPPARAM)
  }
  else
  {
    message("ORA ... ")
    label.prop.gsea(geneset, x, adjacency, threshold = threshold, nperm = nperm,
                    minSize = minSize, maxSize = maxSize, gseaParam = gseaParam,
                    nproc = nproc, BPPARAM = BPPARAM, tmax = tmax, ...)
  }
}
