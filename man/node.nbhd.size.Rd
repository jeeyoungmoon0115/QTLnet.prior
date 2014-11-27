\name{node.nbhd.size}
\alias{node.nbhd.size}

\title{
	Number of moves around a selected node in a network
}
\description{
	Calculates possible number of structures from the current network structure after selecting a node and by performing one of the following moves:  addition [x -> node],  deletion [x -> node],  deleteion [node -> x], or reversion of [x -> node], satisfying \code{max.parents} and acyclicity.
}
\usage{
node.nbhd.size(M, node, max.parents = 3)
}

\arguments{
  \item{M}{
  Adjacency matrix of a phenotype network structure.
}
  \item{node}{
  	A selected node index.
  }
  \item{max.parents}{
  Maximum number of parents.
}
}

\value{
  A list with components
  \item{nbhd.size}{\code{n.deletions + n.additions + n.reversions}}
  \item{n.deletions}{Possible number of structures by a deletion of an edge [x -> node] or [node -> x]}
  \item{n.additions}{Possible number of structures by an addition of an edge [x -> node]}
  \item{n.reversions}{Possible number of structures by a reversion of an edge [x -> node]}
  \item{moves}{A (nbhd.size * 3) data.frame of possible moves. Columns are 'from', 'to', and 'move'. Each row represents the 'move' of an edge ['from' -> 'to'] in the current network.}
}

\note{
 The \code{nbhd.size} component is greater than 0 because  1) if edge addition is not possible, we can always delete an edge [x->node] or [node->x], 2) if there is no edge to delete, in other words, ther is no edge connected with the node, we can always add an edge [x->node] satisfying acyclicity and max.parents.
}

\author{
   Jee Young Moon and Brian S. Yandell
}

\seealso{
\code{\link[qtlnet.prior]{propose.new.node.structure}}	
}

\examples{
mat <- diag(3)
mat[1,2] <- mat[2,3] <- 1
diag(mat) <- 0

node.nbhd.size(mat, 1, 2)
node.nbhd.size(mat, 2, 2)
node.nbhd.size(mat, 3, 2)
	
}

