\name{GO_dme}
\alias{GO_dme}

\title{
Gene Ontology terms with gene ID for Drosophila melanogaster
}
\format{a list of gene sets}

\source{
http://www.go2msig.org/cgi-bin/prebuilt.cgi?taxid=7227
}

\usage{
GO_dme
}

\description{
The data set contains all Gene Ontology terms for Drosophila melanogaster and genes are identified by gene ID. There are 2823 categories.
}


\value{
GO gene sets
}


\author{
Dongmin Jung
}


\examples{
load(system.file("data", "GO_dme.rda", package = "gsean"))
}

