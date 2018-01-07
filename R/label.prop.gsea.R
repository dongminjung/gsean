# label propagation for GSEA
label.prop.gsea <- function(geneset, x, adjacency, threshold = 0.99, nperm = 1000,
                            minSize = 1, maxSize = Inf, gseaParam = 1, nproc = 0,
                            BPPARAM = NULL, ...)
{
  Rownames <- rownames(adjacency)
  Colnames <- colnames(adjacency)
  if(!is.null(Rownames) & !is.null(Colnames))
    if(!identical(Rownames, Colnames))
      stop("the row or column names of the adjacency matrix are not identical")
  if(is.null(Rownames) & is.null(Colnames))
    stop("the row or column names of the adjacency matrix are required")
  
  Names <- Rownames
  if(is.null(Names)) Names <- Colnames
  x <- na.omit(match(x, Names))
  if(length(x) == 0)
    stop("No genes in both the list and the adjacency matrix")
  
  adjacency[adjacency >= threshold] = 1
  adjacency[adjacency < threshold] = 0
  lp <- RANKS::label.prop(adjacency, ind.positives = x, ...)
  scaled.scores <- as.vector(scale(lp$p))
  names(scaled.scores) <- Names
  statistic <- scaled.scores
  
  result.GSEA <- fgsea(geneset, statistic, nperm = nperm, minSize = minSize,
                       maxSize = maxSize, nproc = nproc, gseaParam = gseaParam,
                       BPPARAM = BPPARAM)
  result.GSEA
}
