\name{propose.new.node.structure}
\alias{propose.new.node.structure}

\title{
  Propose a new network with edge addition, deletion or reversion
}
\description{After selecting a node randomly, \code{\link[qtlnet.prior]{node.nbhd.size}} finds possible single moves  by  adding or reversing [x -> node], or by deleting [x -> node] or [node -> x].  This function (\code{\link[qtlnet.prior]{propose.new.node.structure}}) randomly selects a move from the possible moves.  
 }
\usage{
propose.new.structure(M, max.parents = 3, saved.scores, node, neighborM, rev.method = c("node","node.nbhd"), verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
  An adjacency matrix of a phenotype network.
}
  \item{max.parents}{
  Maximum number of parents.
}
  \item{saved.scores}{
Pre-computed scores with all possible parent configuration with \code{max.parents}.
}
  \item{node}{
  	The index of a selected node. 
  }
  \item{neighborM}{
  	Information of possible moves around the selected node. It is a result from \code{\link[qtlnet.prior]{node.nbhd.size}}.
  }
  \item{rev.method}{
  If \code{rev.method="node"}, the move is among the possible single moves from \code{\link[qtlnet.prior]{node.nbhd.size}}. If \code{rev.method="node.nbhd"}, the move for a reversion is enhanced with \code{rev.method='nbhd'} as in Grzegorczyk, M. and Husmeier, D. (2008). 
}
  \item{verbose}{
  	Print out the move.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A list with components
	\item{M}{An adjacency matrix of a new network}
	\item{rev.ratio}{When the move is 'reverse-nbhd', it is a ratio of edge reversion probabilities to the old network  and to the new network.  Otherwise, it is 1.}
	\item{move}{A move among 'add', 'delete', 'reverse', and 'reverse-nbhd'.}
	\item{ne.new}{When the move is 'reverse-nbhd', it is 0 (it is not considered). Otherwise, it is a possible number of structures from the new network by the method \code{rev.method='node'}.}
}

\author{
  Jee Young Moon and Brian S. Yandell
}
\note{
 For \code{rev.method='node'}, after selecting a node randomly with probability \eqn{1/(number of nodes)},  the move is \eqn{1/(neighborM$nbhd.size)}. Depending on the move type, the backward move probability from the new network to the old network is as follows: 1) move was add [x -> node] : backward moves are selecting the node and deleting [node -> x] or selecting x and deleting [node -> x], 2) move was delete [a -> b] where (a,b)=(node,x) or (x,node): selecting node a and adding [a->b], 3) move was reverse [x -> node] : backward move is selecting x and reversing [node -> x].
}

\seealso{
\code{\link[qtlnet.prior]{node.nbhd.size}}	
}

\examples{
mat <- diag(3)
mat[1,2] <- mat[2,3] <- 1
diag(mat) <- 0

node <- 1
max.parents <- 2
neighborM <- node.nbhd.size(mat, node, max.parents)
M.new <- propose.new.node.structure(mat, max.parents=max.parents, saved.scores=NULL, node=node, neighborM=neighborM, rev.method='node', verbose=FALSE)
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
