\name{qtlnet.prior-package}
\alias{qtlnet.prior-package}
\alias{qtlnet.prior}
\docType{package}
\title{
Causal network inference incorporating knowledge
}
\description{
It infers a causal network from gene expression
 and genotypes with prior knowledge	
}
\details{
\tabular{ll}{
Package: \tab qtlnet.prior\cr
Type: \tab Package\cr
Version: \tab 0.6 \cr
Date: \tab 2014-09-16\cr
License: \tab GPL-2 \cr
}
}

\author{
Jee Young Moon and Brian S. Yandell

Maintainer: Jee Young Moon \email{jeeyoung.moon@gmail.com}
}
\references{
  Moon, J.-Y., Chaibub Neto, E., Deng, X.  and Yandell B. S. () Bayesian causal phenotype network
  incorporating genetic variation and biological knowledge. In \emph{Probabilistic Graphical Models Dedicated to Applications in Genetics.} Sinoquet C,
Mourad R, eds. (in review). Manuscript at \url{http://www.stat.wisc.edu/~jmoon/
QTLnet-prior.pdf}
  
  Chiabub Neto, E., Keller, M. P., Attie, A. D. and Yandell, B. S. (2010)
  Causal graphical models in systems genetics: a unified framework
  for joint inference of causal network and genetic archicecture
  for correlated phenotypes.
  \emph{Ann Appl Statist} \bold{4}, 320--339.
  \url{http://dx.doi.org/10.1214/09-AOAS288}

  Grzegorczyk, M. and Husmeier, D. (2008) Improving the structure MCMC sampler
  for Bayesian networks by introducing a new edge reversal move.
  \emph{Mach Learn} \bold{71}, 265--305.
  \url{http://dx.doi.org/10.1007/s10994-008-5057-7}	
}

\seealso{
The main function in the package is \code{\link[qtlnet.prior]{mcmc.qtlnet.prior}}.	
}
\keyword{package}
