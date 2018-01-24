\name{label_prop_gsea}
\alias{label_prop_gsea}

\title{
Over-representaion analysis with the label propagation algorithm
}
\description{
ORA is performed by GSEA with the label propagation algorithm
}
\usage{
label_prop_gsea(geneset, x, adjacency, threshold = 0.99, nperm = 1000,
                minSize = 1, maxSize = Inf, gseaParam = 1, nproc = 0,
                BPPARAM = NULL, ...)
}

\arguments{
  \item{geneset}{
list of gene sets
}
  \item{x}{
set of genes
}
  \item{adjacency}{
adjacency matrix
}
  \item{threshold}{
threshold of correlation for nodes to be considered neighbors (default: 0.99)
}
  \item{nperm}{
number of permutations (default: 1000)
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
fgsea::fgsea, RANKS::label.prop
}

\examples{
data(examplePathways)
data(exampleRanks)
exampleRanks <- exampleRanks[1:100]
geneNames <- names(exampleRanks)
set.seed(1)
x <- sample(geneNames, 10)
adjacency <- diag(length(exampleRanks))
rownames(adjacency) <- geneNames
result.GSEA <- label_prop_gsea(examplePathways, x, adjacency)
}



