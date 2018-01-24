\name{gsean}
\alias{gsean}

\title{
Gene Set Enrichment Analysis with Networks
}
\description{
GSEA or ORA is performed with networks from gene expression data
}
\usage{
gsean(geneset, x, exprs, pseudo = 1, threshold = 0.99, nperm = 1000,
      centrality = function(x) rowSums(abs(x)), weightParam = 1,
      minSize = 1, maxSize = Inf, gseaParam = 1, nproc = 0,
      BPPARAM = NULL, corParam = list(), tmax = 10, ...)
}

\arguments{
  \item{geneset}{
list of gene sets
}
  \item{x}{
Named vector of gene-level statistics for GSEA or set of genes for ORA. Names should be the same as in gene sets.
}
  \item{exprs}{
gene expression data
}
  \item{pseudo}{
pseudo number for log2 transformation (default: 1)
}
  \item{threshold}{
threshold of correlation for nodes to be considered neighbors for ORA (default: 0.99)
}
  \item{nperm}{
number of permutations (default: 1000)
}
  \item{centrality}{
centrality measure, degree centrality or node strength is default
}
  \item{weightParam}{
weight parameter value for the centrality measure, equally weight if weightParam = 0 (default: 1)
}
  \item{minSize}{
minimal size of a gene set (default: 1)
}
  \item{maxSize}{
maximal size of a gene set (default: Inf)
}
  \item{gseaParam}{
GSEA parameter value (default: 1)
}
  \item{nproc}{
see fgsea::fgsea
}
  \item{BPPARAM}{
see fgsea::fgsea
}
  \item{corParam}{
additional parameters for correlation; see WGCNA::cor
}
  \item{tmax}{
maximum number of iterations for label propagtion (default: 10)
}
  \item{...}{
additional parameters for label propagation; see RANKS::label.prop
}
}


\value{
GSEA result
}


\author{
Dongmin Jung
}



\seealso{
exprs2adj, label_prop_gsea, centrality_gsea
}

\examples{
data(examplePathways)
data(exampleRanks)
exampleRanks <- exampleRanks[1:100]
Names <- names(exampleRanks)
exprs <- matrix(rnorm(10*length(exampleRanks)), ncol = 10)
rownames(exprs) <- names(exampleRanks)
set.seed(1)
result.GSEA <- gsean(examplePathways, exampleRanks, exprs)
}



