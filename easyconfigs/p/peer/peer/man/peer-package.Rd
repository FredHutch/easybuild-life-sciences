\name{peer-package}
\alias{peer-package}
\alias{peer}
\alias{PEER}
\alias{PEER_getNc}
\alias{PEER_getPriorAlphaA}
\alias{PEER_getSparsityPrior}
\alias{PEER_init_net_sparse}
\alias{PEER_setPhenoVar}
\alias{PEER_setSigmaOff}
\alias{PEER_calcBound}
\alias{PEER_getNj}
\alias{PEER_getPriorAlphaB}
\alias{PEER_getTolerance}
\alias{PEER_logprob}
\alias{PEER_setPriorAlpha}
\alias{PEER_setSparsityPrior}
\alias{PEER_getAdd_mean}
\alias{PEER_getNk}
\alias{PEER_getPriorEpsA}
\alias{PEER_getVarTolerance}
\alias{PEER_setAdd_mean}
\alias{PEER_setPriorAlphaA}
\alias{PEER_setTolerance}
\alias{PEER_getAlpha}
\alias{PEER_getNmax_iterations}
\alias{PEER_getPriorEpsB}
\alias{PEER_getW}
\alias{PEER_setCovariates}
\alias{PEER_setPriorAlphaB}
\alias{PEER_setPriorEpsB}
\alias{PEER_setVarTolerance}
\alias{PEER_getBounds}
\alias{PEER_getNp}
\alias{PEER_getResidualVars}
\alias{PEER_getX}
\alias{PEER_setNk}
\alias{PEER_setPriorEps}
\alias{PEER_update}
\alias{PEER_getCovariates}
\alias{PEER_getPhenoMean}
\alias{PEER_getResiduals}
\alias{PEER_getZ}
\alias{PEER_setNmax_iterations}
\alias{PEER_setPriorEpsA}
\alias{PEER_getEps}
\alias{PEER_getPhenoVar}
\alias{PEER_getSigmaOff}
\alias{PEER_init_net}
\alias{PEER_setPhenoMean}
\alias{PEER_plotModel}
\alias{PEER_plotQtl}
\docType{package}
\title{
PEER is a collection of Bayesian approaches to infer hidden determinants and their effects from gene expression profiles using factor analysis methods.
~~ package title ~~
}
\description{

PEER is a collection of Bayesian approaches to infer hidden determinants and their effects from gene expression profiles using factor analysis methods. Applications of PEER have detected batch effects and experimental confounders
increased the number of expression QTL findings by threefold allowed inference of intermediate cellular traits, such as transcription factor or pathway activations. 

The PEER model, inference, and applications are described in ``A Bayesian framework to account for complex non-genetic factors in gene expression levels greatly increases power in eQTL studies.'' (PLoS Computational Biology, May 2010)
Joint Genetic Analysis of Gene Expression Data with Inferred Cellular Phenotypes. (PLoS Genetics, January 2011)
and several other projects have successfully used PEER.

This project offers an efficient and versatile C++ implementation of the underlying algorithms with user-friendly interfaces to R and python. To get started using PEER, download the source or binary versions, see the installation instructions, and take a look at the getting started tutorial.

~~ A concise (1-5 lines) description of the package ~~
}
\details{
\tabular{ll}{
Package: \tab peer\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-09-18\cr
License: \tab CC BY-NC-SA 3.0\cr
LazyLoad: \tab yes\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Oliver Stegle, Matias Piipari, Leopold Parts

Maintainer: Matias Piipari <matias.piipari@gmail.com>
~~ The author and/or maintainer of the package ~~
}
\references{

A Bayesian framework to account for complex non-genetic factors in gene expression levels greatly increases power in eQTL studies. (PLoS Computational Biology, May 2010)

Joint Genetic Analysis of Gene Expression Data with Inferred Cellular Phenotypes. (PLoS Genetics, January 2011)

~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
}
\examples{
#TODO: Fill in examples here
}
