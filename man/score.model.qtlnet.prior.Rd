\name{score.model.qtlnet.prior}
\alias{score.model.qtlnet.prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Calculates the BIC score of a network
}
\description{
  It caclulates the BIC score of a network. QTL mapping can be chosen to be done or not.
}
\usage{
score.model.qtlnet.prior(M, saved.scores, cross, addcov, intcov, threshold, nogenotype, verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
   An adjacency matrix of a phenotype network.
}
  \item{saved.scores}{
   Pre-computed BIC score for each child and parent set combination.
}
  \item{cross}{
  Cross object.
}
  \item{addcov}{
  Additive covariates.
}
  \item{intcov}{
  Interactive covariates.
}
  \item{threshold}{
  Threshold for QTL maping.
}
  \item{nogenotype}{
  	Whether QTL mapping is done or not.
  }
  \item{verbose}{
   Verbose.
}
  \item{\dots}{
 Additional arguments.
}
}

\value{
	A list with components
  \item{model.score}{BIC score of the network}
  \item{update.scores}{Updated scores. A list with components 'code','pheno.col', and 'bic'.}
  \item{model.name}{Code name for the network}
}

\examples{
## M
pheno.col=c(1,2,3,4,5)
M=matrix(0, nrow=5, ncol=5)
connected.Y <- c(1,2, 1,3, 1,4,
                 2,5,
                 3,4, 3,5,
                 4,5)
connected.Y <- matrix(connected.Y, ncol=2, byrow=TRUE)

for ( i in seq(nrow(connected.Y))){
  M[connected.Y[i,1], connected.Y[i,2]] <-1
}
	
## Simulate a cross object
cross <- sim.example1.cross(342)
cross <- calc.genoprob(cross, step=0)

threshold <- 3.8 
addcov=NULL; intcov=NULL; max.parents=3;


## Pre-computed score with QTL mapping
saved.scores.QTL <- bic.qtlnet(cross, pheno.col=pheno.col,threshold=threshold, addcov=addcov, intcov=intcov, max.parents=max.parents)
saved.scores.QTL <- bic.join(cross,pheno.col, list(saved.scores.QTL), max.parents=max.parents)

## Pre-computed score without QTL mapping
saved.scores.noQTL <- qtlnet.prior:::bic.qtlnet.pheno(cross, pheno.col=pheno.col,threshold=threshold, addcov=addcov, intcov=intcov, max.parents=max.parents)
saved.scores.noQTL <- bic.join(cross,pheno.col, list(saved.scores.noQTL), max.parents=max.parents)

## Score of a network      	
score.model.qtlnet.prior(M, saved.scores=saved.scores.QTL, cross=cross, addcov=NULL, intcov=NULL, threshold=3.8, nogenotype=FALSE, verbose=TRUE)
score.model.qtlnet.prior(M, saved.scores=saved.scores.noQTL, cross=cross, addcov=NULL, intcov=NULL, threshold=3.8, nogenotype=TRUE, verbose=TRUE)
	
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
