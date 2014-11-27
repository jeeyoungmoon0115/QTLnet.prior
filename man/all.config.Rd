\name{all.config}
\alias{all.config}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Columnwise calculation of |B-G| for each parent configuration
}
\description{
|B-G|  is calculated for each column i: |B[,i] - G[,i]|. For each child node/column i, all possible parent configurations up to \code{max.parents.partition} are considered and |B[,i] - G[,i]| is calculated for each parent configuration k.
}
\usage{
all.config(n.pheno, max.parents.partition, B)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.pheno}{
  	Number of phenotypes.
}
  \item{max.parents.partition}{
   Maximum number of parents to consider for partition function calculation. 
}
  \item{B}{
  	A list of biological knowledge matrices.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A list with components
	\item{e.value}{It is an array of size (n.pheno, count.parents.config, number of types of knowledge). \code{e.value[i,k,j]} is \eqn{| B[,i] - G[,i] | } for child node i, parent configuration k, and j-th knowledge. It is calcualted by \eqn{sum(1-B[[j]][k,i]) + sum(B[[j]][-k,i])}.}
	\item{count.parents.config}{Total number of parent configurations from no parent up to \code{max.parents.parition} parents for a node.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
  \code{\link[qtlnet.prior]{partition.eff}}
}
\examples{
B <- runif(25, 0, 1)
B <- matrix(B, ncol=5)
diag(B) <- 0
B <- list(B)

ex.config <- all.config(n.pheno=5, max.parents.partition=3, B)  
}
