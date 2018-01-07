# expression to adjacency matrix
exprs2adj <- function(x, pseudo = 1, ...)
{
  index <- which(rowSums(x) == 0)
  if(length(index) > 0)
  {
    y <- x[-index,]
  }
  else
  {
    y <- x
  }
  
  # log2 transform from GEO2R 
  qx <- as.numeric(quantile(y, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=TRUE))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC && min(y, na.rm = TRUE) >= 0)
  {
    y <- log2(y + pseudo)
    cat("log2 transformed for correlation")
  }
  A <- WGCNA::cor(t(y), ...)
  A
}
