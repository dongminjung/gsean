# centrality measure for GSEA
centrality_gsea <- function(geneset, x, adjacency, pseudo = 1, nperm = 1000,
                            centrality = function(x) rowSums(abs(x)),
                            weightParam = 1, minSize = 1, maxSize = Inf,
                            gseaParam = 1, nproc = 0, BPPARAM = NULL)
{
  Rownames <- rownames(adjacency)
  Colnames <- colnames(adjacency)
  if(!is.null(Rownames) & !is.null(Colnames))
    if(!identical(Rownames, Colnames))
      stop("the row or column names of the adjacency matrix are not identical")
  if(is.null(Rownames) & is.null(Colnames))
    stop("the row or column names of the adjacency matrix are required")
  
  x.names <- names(x)
  Names <- Rownames
  if(is.null(Names)) Names <- Colnames
  
  overlap <- na.omit(intersect(x.names, Names))
  if(length(overlap) == 0)
    stop("No genes in both the names of statistic and the adjacency matrix")
  
  x <- x[overlap]
  if(is.null(Rownames))
    rownames(adjacency) <- Names
  if(is.null(Colnames))
    colnames(adjacency) <- Names
  adjacency <- adjacency[overlap, overlap]
  
  # log2 transform from GEO2R 
  qx <- as.numeric(quantile(x, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=TRUE))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC && min(x, na.rm = TRUE) >= 0)
  {
    x <- log2(x + pseudo)
    cat("log2 transformed for statistic")
  }
  
  statistic <- x*(centrality(adjacency)^weightParam)
  result.GSEA <- fgsea(geneset, statistic, nperm = nperm, minSize = minSize,
                       maxSize = maxSize, nproc = nproc, gseaParam = gseaParam,
                       BPPARAM = BPPARAM)
  result.GSEA
}
