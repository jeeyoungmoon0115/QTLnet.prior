\name{partition.eff}
\alias{partition.eff}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Calculates the log of partition function
}
\description{
	It calculates the log of partition function in an approximate way (normalizing term in prior probability of phenotype networks) specified by beta.   The parition function is defined to be \eqn{\sum_{G_Y} exp(- \sum_j beta_j * |B[[j]] - G_Y|).}  Instead, for each child node i, it first calculates the \eqn{\sum_j |B[[j]][,i] - G[,i]|}  over all parent configurations: \eqn{\sum_{k in node i's parent configuration} exp(- \sum_j beta_j * [(1-B[[j]][k,i]) + (B[[j]][-k,i])] ).}  Then, the product of the values across nodes are used as the parition function. This function \code{partition.eff} returns the log of the calculated partition function.
}
\usage{
partition.eff(beta, n.pheno, B, count.parents.config, e.value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
  A vector of beta values. beta is the same length of list \code{B}.
}
  \item{n.pheno}{
  Number of phenotypes.
}
  \item{B}{
    A list of biological knowledge matrices.
}
  \item{count.parents.config}{
  Total number of possible parents' configurations up to \code{max.parents.partition} parents for a node. See \code{\link[qtlnet.prior]{all.config}}.
}
  \item{e.value}{
  	An array of (n.pheno, count.parents.config, number of types of knowledge). \code{e.value[i,k,j]} is |B[[j]][,i] - G_Y[,i]| with G_Y having the k-th parent configuration for the child node i, and j-th knowledge.  See \code{\link[qtlnet.priro]{all.config}}.
%%     ~~Describe \code{e.value} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A numeric value, log of the partition function.
}
\references{
    Werhli, A. V. and Husmeier, D. (2007) Reconstructing gene regulatory networks 
  with Bayesian networks by combining expression data with multiple sources of
  prior knowledge. \emph{Stat. Appl. Genet. Molec.} \bold{6}:15.
}
\author{
  Jee Young Moon and Brian S. Yandell
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link[qtlnet.prior]{all.config}}, 
  \code{\link[qtlnet.prior]{score.prior}}
}
\examples{
B <- runif(25, 0, 1)
B <- matrix(B, ncol=5)
diag(B) <- 0
B <- list(B)

ex.config <- all.config(n.pheno=5, max.parents.partition=3, B)  
partition.eff(beta=1.2, n.pheno=5, B, ex.config[['count.parents.config']], ex.config[['e.value']])
}

