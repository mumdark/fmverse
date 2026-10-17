% Generated documentation for the fmverse Geneformer registry.
\name{fmverse_geneformer_checkpoints}
\alias{fmverse_geneformer_checkpoints}
\title{List Supported Geneformer Checkpoints}
\usage{
fmverse_geneformer_checkpoints()
}
\value{
A data frame describing the four current official Geneformer core checkpoints,
including model family, parameter count, training corpus, input size, dictionary
release, and whether special boundary tokens are used.
}
\description{
Returns the version-specific checkpoint and preprocessing specifications used
by the Geneformer adapter.
}
\examples{
fmverse_geneformer_checkpoints()
}
