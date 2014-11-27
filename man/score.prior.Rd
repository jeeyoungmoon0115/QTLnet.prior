\name{score.prior}
\alias{score.prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Log of the numerator in prior probability of phenotype networks
}
\description{
	The prior probablity of phenotype networks is defined to be \eqn{P(G_Y | B, beta) \propto exp(-\sum_j beta[j] * |B[[j]]-M|)}. The normalizing term (denominator) for this probability is called a partition function. This function \code{score.prior} calculates the log of the numerator and \code{\link[qtlnet.prior]{partition.eff}} approximately calculates the log of the denominator.
}
\usage{
score.prior(M, beta, B)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
  Adjacency matrix.
}
  \item{beta}{
   A vector of beta values. \code{beta[j]} corresponds to \code{B[[j]]}.
}
  \item{B}{
  A list of matrices to encode biologicla knowledge. Each matrix (\code{B[[j]]}) corresponds to one type of biological knowledge.
}
}
\details{
If beta is close to 0, the prior is close to a uniform distribution and hence the contribution of biological knowledge is negligible. If \eqn{beta -> infinity}, the prior puts the most of probability on the structure closest to B.   
}
\value{
	It returns a numeric value : \eqn{-1 * \sum_j beta[j] * |M - B[[j]]|}. 

}

\author{
 Jee Young Moon and Brian S. Yandell
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\references{
	  Werhli, A. V. and Husmeier, D. (2007) Reconstructing gene regulatory networks 
  with Bayesian networks by combining expression data with multiple sources of
  prior knowledge. \emph{Stat. Appl. Genet. Molec.} \bold{6}:15.
}

\seealso{
\code{\link[qtlnet.prior]{partition.eff}}
}
\examples{
M.new <- diag(5)
M.new[1,2] <- M.new[2,3] <- 1
diag(M.new) <- 1

beta <- 1.2

B <- M.new
B[1,2] <- 0.8
B[2,3] <- 0.7
B <- list(B)

M <- diag(5)
M[3,4] <- M[4,5] <- 1
diag(M) <- 0
	
score.prior(M.new, beta, B)	
score.prior(M, beta, B)	

}
