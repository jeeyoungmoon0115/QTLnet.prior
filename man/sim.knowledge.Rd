\name{sim.knowledge}
\alias{sim.knowledge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Simulate a knowledge matrix
}
\description{
  Simulate a knowledge matrix
}
\usage{
sim.knowledge(B.N = 1, delta = 0.5, A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{B.N}{
  Number of types of biological knowledge
}
  \item{delta}{
  A vector of accuracy for each knowledge type. Delta value is between -0.5 and 0.5. The higher the delta is, the more accurate the knowledge will be. 
}
  \item{A}{
  True adjacency matrix of the network 
}
}
\details{
	Each element correponding to 1 in the true adjacency matrix \code{A} is simulated from a truncated normal distribution between 0 and 1 with mean \eqn{0.5 + delta} and sigma=0.1. Each element correponding to 0 in the true adjacency matrix \code{A} is simulated from a truncated normal distribution between 0 and 1 with mean \eqn{0.5 - delta} and sigma=0.1. 
}
\value{
	A list with B.N components and each component is a matrix of the same size as \code{A}. Each component corresponds to one type of biological knowledge.
}

\author{
Jee Young Moon and Brian S. Yandell
}

\examples{
A <- diag(3)
A[1,2] <- A[2,3] <- 1
diag(A) <- 0
	
sim.knowledge(B.N=2, delta=0.3, A)
sim.knowledge(B.N=2, delta=c(0.3, -0.3), A)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
