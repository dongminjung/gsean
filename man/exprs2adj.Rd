\name{exprs2adj}
\alias{exprs2adj}

\title{
Convert gene expression data to adjacency matrix by using correlation coefficients
}
\description{
A biological network is constructed from gene expression data and it is used for Gene Set Enrichment Analysis.
}
\usage{
exprs2adj(x, pseudo = 1, ...)
}

\arguments{
  \item{x}{
gene expression data
}
  \item{pseudo}{
pseudo number for log2 transformation (default: 1)
}
  \item{...}{
additional parameters for correlation; see WGCNA::cor
}
}


\value{
adjacency matrix
}


\author{
Dongmin Jung
}



\seealso{
fgsea::fgsea, WGCNA::cor
}

\examples{
data(exampleRanks)
Names <- names(exampleRanks)
exprs <- matrix(rnorm(10*length(exampleRanks)), ncol = 10)
adjacency <- exprs2adj(exprs)
}



