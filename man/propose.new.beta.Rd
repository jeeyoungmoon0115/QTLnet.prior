\name{propose.new.beta}
\alias{propose.new.beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Propose a new beta
}
\description{
  Propose a new beta. First, it samples \code{beta.new} from \code{U[beta.old-1, beta.old+1]}. Second, if \code{beta.new} is negative, use \code{-1 * beta.new}. This proposal satisfies the following ratio: \eqn{P(beta.new | beta.old) / P(beta.old | beta.new) =1}.}
  
\usage{
propose.new.beta(beta.old, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta.old}{
  	Old beta
}
  \item{lambda}{
  lambda is not used.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A positive value for a new beta.
}
\author{
    Jee Young Moon and Brian S. Yandell
}

\examples{
propose.new.beta(0.5,1)
}
