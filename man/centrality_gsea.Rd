\name{centrality_gsea}
\alias{centrality_gsea}

\title{
Gene Set Enrichment Analysis with centrality measure
}
\description{
GSEA is performed with centrality measure
}
\usage{
centrality_gsea(geneset, x, adjacency, pseudo = 1, nperm = 1000,
                centrality = function(x) rowSums(abs(x)),
                weightParam = 1, minSize = 1, maxSize = Inf,
                gseaParam = 1, nproc = 0, BPPARAM = NULL)
}

\arguments{
  \item{geneset}{
list of gene sets
}
  \item{x}{
Named vector of gene-level statistics. Names should be the same as in gene sets.
}
  \item{adjacency}{
adjacency matrix
}
  \item{pseudo}{
pseudo number for log2 transformation (default: 1)
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
}


\value{
GSEA result
}


\author{
Dongmin Jung
}



\seealso{
fgsea::fgsea
}

\examples{
data(examplePathways)
data(exampleRanks)
exampleRanks <- exampleRanks[1:100]
adjacency <- diag(length(exampleRanks))
rownames(adjacency) <- names(exampleRanks)
set.seed(1)
result.GSEA <- centrality_gsea(examplePathways, exampleRanks, adjacency)
}



