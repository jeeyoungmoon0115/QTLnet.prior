\name{score.beta}
\alias{score.beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Log of the prior probability of beta
  }
\description{
	The hyperprior probablity of beta is defined to be \eqn{P(beta | lambda) = lambda e^(-1 * lambda * beta).}{P(beta | lambda) = lambda exp(-1 * lambda * beta),} where \eqn{lambda} is a fixed value. This function \code{score.beta} returns \eqn{-lambda * beta} without a constant term as a score for the prior of beta. 
}
\usage{
score.beta(beta, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
  Beta.
}
  \item{lambda}{
   lambda.
}

}

\value{
	It returns a numeric value : \eqn{-1 * lambda * beta}. 

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
 Jee Young Moon and Brian S. Yandell
}
\note{
%%  ~~further notes~~
}


