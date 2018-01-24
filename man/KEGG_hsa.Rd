\name{KEGG_hsa}
\alias{KEGG_hsa}

\title{
KEGG pathways with gene symbol for human
}
\format{a list of gene sets}

\source{
http://software.broadinstitute.org/gsea/msigdb/collections.jsp
}

\usage{
KEGG_hsa
}

\description{
The data set contains 186 KEGG pathways for Drosophila melanogaster and genes are identified by gene symbol.
}


\value{
KEGG gene sets
}


\author{
Dongmin Jung
}


\examples{
load(system.file("data", "KEGG_hsa.rda", package = "gsean"))
}


