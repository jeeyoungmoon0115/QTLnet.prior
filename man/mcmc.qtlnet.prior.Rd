\name{mcmc.qtlnet.prior}
\alias{mcmc.qtlnet.prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 MCMC run for causal network inference with a prior knowledge
}
\description{
	It runs an MCMC to infer a causal network of genotypes and phenotypes 
	incorporating knowledge as a prior.
}
\usage{
mcmc.qtlnet.prior(cross, B = NULL, nogenotype=FALSE, pheno.col=NULL, 
  addcov = NULL, intcov = NULL, 
  nSamples = 3000, thinning = 10, burnin = 0.1, random.seed = NULL,
  max.parents = 3, M0 = NULL, init.edges = 0, beta0 = NULL, lambda = NULL,
  threshold=NULL, n.perm=1000, alpha=0.05, step=0, method = "hk", 
  saved.scores = NULL, rev.method = c("node.nbhd","node","nbhd", "node.edge", "single"), 
  verbose = FALSE, max.parents.partition = 3,  filename = "tmp.RData", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{
   An object of class \code{cross}. It contains genotypes and phenotypes. See \code{\link[qtl]{read.cross}}.
}
  \item{B}{
   A matrix if only one type of biological knowledge is incorporated or a list of matrices for biological knowledge. If \code{B = NULL}, biological knowledge is not used in the causal network inference. Otherwise, each matrix in the list encodes one type of biological knowledge. Each matrix is a square matrix where (i,j)-th element is the probability of edge presence (i -> j) supported by the corresponding type of biological knowledge. Elements should be between 0 and 1. 
}
  \item{nogenotype}{
If genotype information is used, \code{nogenotype = FALSE}. If genotype information is not used, 	 \code{nogenotype = TRUE}.  	
  } 
  \item{pheno.col}{
  	  Phenotype identifiers from \code{cross} object to be included in the causal network inference. May be numeric, logical or character. If \code{pheno.col=NULL}, all the phenotypes in the cross object will be used.
}


  \item{addcov}{
   Additive covariates for each phenotype. If there is no additive covariate,  \code{addcov = NULL}.
}
  \item{intcov}{
   Interactive covariates for each phenotype. If there is no interactive covariate,  \code{intcov = NULL}.
}


  \item{nSamples}{
   Number of MCMC samples to record.
}
  \item{thinning}{
   Every \code{thinning} is recorded after \code{burnin} period.
}
  \item{burnin}{
   Ratio of initial MCMC iterations to drop.
}
  \item{random.seed}{
   Random seed.
}


  \item{max.parents}{
  Maximum number of parents allowed in the network structure. 
}
  \item{M0}{
   Initial phenotype network structure satisfying \code{init.edges} and \code{max.parents} in the MCMC run. It is a \sQuote{number of phenotypes} by \sQuote{number of phenotypes} matrix.
}
  \item{init.edges}{
  	Number of edges in the initial network.
}
  \item{beta0}{
   Initial beta value in the MCMC run. \code{beta} is a hyperparameter in the prior distribution for phenotype networks which controls the contribution of biological knowledge. The prior distribution for phenotype networks is defined to be \eqn{P(G_Y | B, beta) \propto e^(- beta |G_Y - B|).}{P(G_Y | B, beta) \propto exp(- beta |G_Y - B|).}
}
  \item{lambda}{
  	The scale parameter in the hyperprior distribution for beta.  It sets the hyperprior for beta to be  \eqn{P(beta | lambda) = lambda e^(-lambda * beta).}{P(beta | lambda) = lambda exp(-lambda * beta).}
}


  \item{threshold}{
  	 A scalar or a vector of thresholds for QTL identification. If it is a scalar, it will be replicated as a vector of  \code{length(pheno.col)}. Otherwise, it should be the same length as the number of phenotypes. The same threshold is used for autosomes and X-chrosomome for this version. If genotypes are not included in network inference (\code{nogenotype=TRUE}), \code{threshold} will be a vector of \code{NA}.
}

  \item{n.perm}{
  	If \code{threshold} is \code{NULL}, a permutation test to get a vector of thresholds for QTL identification is performed. \code{n.perm} is the number of permutations to perform to get an \code{alpha}-level threshold for QTL detection for each phenotype. If there is \code{addcov} or \code{intcov}, the permutation test takes into account of it.
  }
  \item{alpha}{
  	Significance cut-off for thresholds for QTL identification. 
  }
  \item{step}{
  	Genotype probabilities are calculated at the maximum \code{step} (cM) between pseudomarkers. When \code{step = 0}, genotype probabilities are calculated only at the marker positions.
  }
    \item{method}{
  	A method used for QTL interval mapping. A default is Haley-Knott regression. See \code{\link[qtl]{scanone}}.
 }



  \item{saved.scores}{
  	Pre-computed scores with all possible parent configurations upt to \code{max.parents}. If it is \code{NULL}, it will be calculated inside the function \code{mcmc.qtlnet.prior}.
}
  \item{rev.method}{
  	A proposal method for a new network structure.  If \code{rev.method = 'single'}, a proposal is among 'add', 'delete', and 'reverse' of an edge. If \code{rev.method = 'nbhd'}, it is an updated version of \code{rev.method='single'} in reversing an edge as described in Grzegorczyk and Husmeier (2008).
  If \code{rev.method='node.edge'}, it is an extension of the method in Grzegorczyk and Husmeier (2008). In proportional to the score between a parent set and a node, it randomly selects 1) a node and its new parent set, 2) an existing edge and new parent sets for its two nodes, without making an edge between them, or 3)  an existing edge and new parent sets for its two nodes, reversing the edge between them.  
  If \code{rev.method='node'}, it randomly selects a node and moves to the new structure by adding or reversing [x -> node], or deleting [x -> node] or [node -> x].  If \code{rev.method='node.nbhd'}, it is an mix of 'node' and 'nbhd' methods.
  	}


  \item{verbose}{
  	If \code{TRUE} or \code{1}, print iteration numbers and permutations. If \code{verbose=2} or {4}, print MCMC moves. If \code{verbose=3} or {4}, plot BIC. 
}
  \item{max.parents.partition}{
  	Maximum number of parents to be used to calculate the partition function (normlizing term in the prior probability of phenotype network structures).
}
  \item{filename}{
  The name of a file to save the output.
}
  \item{\dots}{
  	Additional arguments.
}
}
\details{
It runs an MCMC to infer a causal network of genotypes and phenotypes, incorporating biological knowledge \code{B}. 

Arguments related with using data: \code{nogenotype}, \code{B}. This function can run with four conditions as described in the example below -- with or without genotypes, with or without biological knowledge.

Arguments related with initial points in MCMC: \code{M0}, \code{init.edges}, \code{beta0}.

Arguments related with MCMC running length: \code{nSamples}, \code{thinning}, \code{burnin}.

Arguments related with QTL mapping when using genotypes: \code{threshold}, \code{n.perm}, \code{alpha}, \code{step}, \code{method}.

Arguments to constrain the network structures due to the super-exponential size of network structure space: \code{max.parents}.

Arguments for network structure proposal: \code{rev.method}. There are five ways to propose a new network structure. 
 }
\value{
	List of class \code{qtlnet} with additional values. See \code{\link[qtlnet]{mcmc.qtlnet}}. Additional values are
	\item{mav}{Model average of \code{M} for \code{nSamples} of MCMC samples}
	\item{cont.accept}{Count of acceptance or rejection for new structure proposals (addition, deletion, reverse)}
	\item{post.beta}{Posterior beta}
	\item{post.partition}{Partition function calculation (normalizing term in the prior probability of phenotype networks) specified by beta value} 
	\item{beta.freq.accept}{Acceptance frequency of beta}
	\item{post.prob}{Proportional to log of posterior probability}
	
	
	Additional attributes from input
	\item{B}{A list of matrices for biological knowledge}
	\item{nogenotype}{Indicator whether genotypes are used}
	\item{lambda}{The scale parameter in the hyperprior distribution for beta.}
}
\references{
  Moon J.-Y., Chaibub Neto, E., Deng, X. and Yandell B. S. () Bayesian causal phenotype network
  incorporating genetic variation and biological knowledge. In \emph{Probabilistic Graphical Models Dedicated to Applications in Genetics.} Sinoquet C,
Mourad R, eds. (in review). Manuscript at \url{http://www.stat.wisc.edu/~jmoon/
QTLnet-prior.pdf}
  
  Chiabub Neto, E., Keller, M. P., Attie, A. D. and Yandell, B. S. (2010)
  Causal graphical models in systems genetics: a unified framework
  for joint inference of causal network and genetic archicecture
  for correlated phenotypes.
  \emph{Ann Appl Statist} \bold{4}, 320--339.
  \url{http://dx.doi.org/10.1214/09-AOAS288}

  Husmeier, D. (2003) Sensitivity and specificity of inferring genetic regulatory interactions
  from microarray experiments with dynamic Bayesian networks. \emph{Bioinformatics} \bold{19}:2271--2282
  
  Werhli, A. V. and Husmeier, D. (2007) Reconstructing gene regulatory networks 
  with Bayesian networks by combining expression data with multiple sources of
  prior knowledge. \emph{Stat. Appl. Genet. Molec.} \bold{6}:15.

  Grzegorczyk, M. and Husmeier, D. (2008) Improving the structure MCMC sampler
  for Bayesian networks by introducing a new edge reversal move.
  \emph{Mach Learn} \bold{71}, 265--305.
  \url{http://dx.doi.org/10.1007/s10994-008-5057-7}	

}
\author{
  Jee Young Moon and Brian S. Yandell
  
  Maintainer: Jee Young Moon \email{jeeyoung.moon@gmail.com}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For the class of \code{qtlnet}:	\code{\link[qtlnet]{mcmc.qtlnet}}. 

For the \code{cross} object and QTL mapping: \code{\link[qtl]{read.cross}}, \code{\link[qtl]{scanone}}.

For the prior probability of a phenotype network: \code{\link[qtlnet.prior]{score.prior}}, \code{\link[qtlnet.prior]{partition.eff}}, \code{\link[qtlnet.prior]{all.config}}. 

For the prior probabiility of beta: \code{\link[qtlnet.prior]{score.beta}}.

For the BIC score of a network: \code{\link[qtlnet.prior]{score.model.qtlnet.prior}}.

For the proposal of new network structures by \code{rev.method='node.nbhd'}: \code{\link[qtlnet.prior]{propose.new.node.structure}}, \code{\link[qtlnet.prior]{node.nbhd.size}}.

For the generation of simulated knowledge: \code{\link[qtlnet.prior]{sim.knowledge}}.

%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## A network of 5 phenotypes with following causal relations:
##    1->2, 1->3, 1->4, 2->5, 3->4, 3->5, 4->5.
## i->j is denoted by A[i,j]=1; otherwise, A[i,j]=0
pheno.col=c(1,2,3,4,5)
A=matrix(0, nrow=5, ncol=5)
connected.Y <- c(1,2, 1,3, 1,4,
                 2,5,
                 3,4, 3,5,
                 4,5)
connected.Y <- matrix(connected.Y, ncol=2, byrow=TRUE)

for ( i in seq(nrow(connected.Y))){
  A[connected.Y[i,1], connected.Y[i,2]] <-1
}

## Number of biological knowledge types
B.N <- 1
## Accuracy of biological knowledge
delta <- c(0.2)
## scanone threshold for QTL detection
threshold = 3.83

## MCMC simulation parameters
# final number of MCMC samples
nSamples = 1000
# MCMC samples will  
thinning = 10
# Maximum number of parents
max.parents = 3

## Simulate a cross object
cross <- sim.example1.cross(342)
cross <- calc.genoprob(cross, step=0)
## Simulate a list of knowledge matrices
B <- sim.knowledge(B.N, delta, A)

## Save scanone result conditonal on the possible parental configurations
##    under the max.parents
saved.scores<-bic.qtlnet(cross, pheno.col=pheno.col,threshold=threshold, max.parents=max.parents)
saved.scores<-bic.join(cross,pheno.col, list(saved.scores), max.parents=max.parents)

## Case 1: Run MCMC to infer a causal network of genotypes and phenotypes 
## incorporating biological knowledge
outBG <- mcmc.qtlnet.prior(cross, B=B, nogenotype=FALSE, pheno.col=pheno.col, 
              nSamples=nSamples, thinning=thinning, random.seed=243, 
              init.edges = 0, threshold=threshold, saved.scores=saved.scores,
              rev.method='nbhd')
              
## Case 2: Run MCMC to infer a causal network of genotypes and phenotypes
## (no biological knowledge) :  Chaibub Neto, E. et al. (2010)
outG <- mcmc.qtlnet.prior(cross, B=NULL, nogenotype=FALSE, pheno.col=pheno.col, 
              nSamples=nSamples, thinning=thinning, random.seed=472, 
              init.edges = 0, threshold=threshold, saved.scores=saved.scores,
              rev.method='nbhd')

## Case 3: Run MCMC to infer a causal network of phenotypes 
## incorporating biological knowledge (no genotypes) :  Werhli, A. V. and Husmeier, D. (2007)
outB <- mcmc.qtlnet.prior(cross, B=B, nogenotype=TRUE, pheno.col=pheno.col, 
              nSamples=nSamples, thinning=thinning, random.seed=6821, 
              init.edges = 0, threshold=NULL, saved.scores=NULL,
              rev.method='nbhd')

## Case 4: Run MCMC to infer a causal network of phenotypes 
## (no genotypes, no biological knowledge)
out <- mcmc.qtlnet.prior(cross, B=NULL, nogenotype=TRUE, pheno.col=pheno.col, 
              nSamples=nSamples, thinning=thinning, random.seed=576, 
              init.edges = 0, threshold=NULL, saved.scores=NULL,
              rev.method='nbhd')              

}

\keyword{qtlnet.prior}
