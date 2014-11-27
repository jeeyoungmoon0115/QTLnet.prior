#######################################################################
## 2014/07/17
##  Copyright (C) 2014 Jee Young Moon and Brian S. Yandell
## This code is under the GNU.
##
## Note: mcmc.qtlnet.prior() is modified from qtlnet::mcmc.qtlnet()
## to set a prior probability from biological knowledge and
## a hyperparameter for a tuning parameter.
##
## Note: nbhd.size() is slighly modified from qtlnet:::nbhd.size()
##
## Functions: mcmc.qtlnet.prior, nbhd.size, check.additions, check.reversions.edge, check.cycle,
##   propose.new.beta, score.prior, score.beta, all.config,
##  partition.eff
######################################################################
mcmc.qtlnet.prior <- function(cross, B = NULL, nogenotype=FALSE, pheno.col=NULL,
                        addcov=NULL, intcov=NULL,
                        nSamples = 3000, thinning=10, burnin=0.1, random.seed=NULL,
                        max.parents = 3,  M0 = NULL, init.edges = 0, beta0=NULL, lambda=NULL,
                        threshold=NULL, n.perm=1000, alpha=0.05, step=0, method="hk",
                        saved.scores = NULL, rev.method = c("node.nbhd","node","nbhd", "node.edge", "single"),
                        verbose = FALSE, max.parents.partition=3, filename='tmp.RData',...)
{
	
  ## Verbose: 1 or TRUE: saved count; 2: MCMC moves; 3: plot BIC; 4: 2&3.
  
  ## Check input parameters.
  rev.method <- match.arg(rev.method)

  ## Random number generator seed.
  if(!is.null(random.seed)) {
    if(!is.numeric(random.seed))
      stop("random seed must be numeric")
    set.seed(random.seed) 
  }

  ## Burnin must be between 0 and 1.
  if(is.logical(burnin))
    burnin <- ifelse(burnin, 0.1, 0)
  if(burnin < 0 | burnin > 1)
    stop("burnin must be between 0 and 1")

  ## Adjust phenotypes and covariates to be numeric.
  if (is.null(pheno.col))  pheno.col <- seq(nphe(cross))
  cross <- qtlnet:::adjust.pheno(cross, pheno.col, addcov, intcov)
  pheno.col <- cross$pheno.col
  pheno.names <- cross$pheno.names
  addcov <- cross$addcov
  intcov <- cross$intcov
  cross <- cross$cross

  n.pheno <- length(pheno.col)

  ## Calculate genotype probabilities if not already done.
  if (!("prob" %in% names(cross$geno[[1]]))) {
    #warning("First running calc.genoprob.")
    cross <- calc.genoprob(cross, step=step)
  }

  ## LOD threshold by phenotype.
  max.parents <- min(max.parents, n.pheno-1)
  if(!nogenotype){
    ## when genotypes are included in the network
    if (is.null(threshold)){
      cat("## Thresholds for QTL identification are not provided.\n## Calculating ",alpha," threshold for each phenotype.\n")
      permo <- scanone(cross, pheno.col=pheno.col, method=method, addcovar=addcov, intcovar=intcov, n.perm=n.perm, perm.Xsp=FALSE, verbose=verbose)
      threshold <- c(summary(permo, alpha=alpha))
    }
  
    if(length(threshold) == 1)
    threshold <- rep(threshold, n.pheno)
    if(length(threshold) != n.pheno) stop("threshold must have same length as pheno.col")
  }
  else{
      ## when genotypes are not included in the network
      threshold <- rep(NA, n.pheno)
  }
  

  ## Initial network matrix M0.
  if(is.null(M0))  M0 <- init.qtlnet(pheno.col, max.parents, init.edges)
  if(nrow(M0) != n.pheno | ncol(M0) != n.pheno)
    stop("M0 must be square matrix the size of pheno.col")

  ## Check B.
  if (class(B)=='matrix') B <- list(B)
  if (class(B)=='list'){
      for(i in seq(length(B))){
          if (ncol(B[[i]]) != n.pheno) stop("number of columns in B[[i]] should be n.pheno")
          if (nrow(B[[i]]) != n.pheno) stop("number of rows in B[[i]] should be n.pheno")
          if (any(!is.numeric(B[[i]]))) stop("elements in B[[i]] should be numeric.")
          if (any(diag(B[[i]]) !=0)) stop("diagonal in B[[i]] should be 0.")
          if (any(B[[i]]>1 | B[[i]]<0)) stop("elemtns in B[[i]] should be between 0 and 1.")
      }
  }
  
  ## Initial beta
  if (!is.null(B))  {
    if (is.null(lambda))  lambda=rep(1,length(B)) #lambda = energy.diff(B)
    else if (length(lambda)!=length(B)) lambda=rep(lambda[1],  length(B))
  }
  else lambda=NULL
  if(is.null(beta0)) beta0<- rexp(length(B), lambda)


  ## saved.scores
  if(is.null(saved.scores)){
    #rev.method <- "single"
    if (!nogenotype){
      ## when genotypes are included in the network
      saved.scores<-bic.qtlnet(cross, pheno.col=pheno.col,threshold=threshold, addcov=addcov, intcov=intcov, max.parents=max.parents)
      saved.scores<-bic.join(cross,pheno.col, list(saved.scores), max.parents=max.parents)
    }
    else{
      ## when genotypes are not included in the network
      saved.scores<-bic.qtlnet.pheno(cross, pheno.col=pheno.col,threshold=threshold, addcov=addcov, intcov=intcov, max.parents=max.parents)
      saved.scores<-bic.join(cross,pheno.col, list(saved.scores), max.parents=max.parents)
    }
  }
  saved.scores <- qtlnet:::make.saved.scores(pheno.names, max.parents,
       saved.scores = saved.scores, verbose = verbose, ...)


  ## Create output objects to be used for recording
  mav <- Mav <- matrix(0, n.pheno, n.pheno)
  n.burnin <- burnin * nSamples * thinning
  
  post.bic <- rep(NA,nSamples)
  post.model <- rep(NA,nSamples)
  all.bic <- rep(NA,nSamples)
  post.beta <- matrix(NA,nrow=nSamples,ncol=length(B))
  post.partition <- rep(NA, nSamples)
  post.prob <- rep(NA, nSamples)
  #B.mean<-lapply(B, function(x) sum(x)/n.pheno)
  
  
  M.old <- M0
  ne.old <- qtlnet:::nbhd.size(M.old, max.parents)[[1]]

  aux.new <- score.model.qtlnet.prior(M.old, saved.scores, cross, addcov, intcov,
                         threshold, nogenotype=nogenotype, verbose, method = method, ...)
  
  bic.old <- aux.new$model.score
  tmp <- aux.new$update.scores
  codes <- dimnames(saved.scores)[[1]]
  if(!is.null(tmp$bic)) {
    index <- nrow(saved.scores)
    index <- match(as.character(tmp$code), codes) + index * (tmp$pheno.col - 1)
    saved.scores[index] <- tmp$bic
  }
  model.old <- aux.new$model.name


  beta.new <- beta.old <-beta0
  Mprior.old <- score.prior(M.old, beta.old, B)
  prior.beta.old <- score.beta(beta.old, lambda)
  if(is.null(B)) partition.old <- NA
  else {
    pa.config<-all.config(n.pheno, max.parents.partition, B)
    count.parents.config <- pa.config[[2]]
    e.value <- pa.config[[1]]
    partition.old <- partition.eff(beta.old, n.pheno, B, count.parents.config, e.value)
    #partition.old <- partition(beta.old, n.pheno, B, max.parents.partition, B.mean)
  }

  if (!is.null(B)) post.prob.old <- bic.old - 2*Mprior.old + 2 * partition.old


  if(verbose) {
    cat("\n")
    if(verbose > 2)
      plot(c(1,nSamples), bic.old * c(0.5,1.2), type = "n",
           xlab = "sample", ylab = "BIC")
  }


  k <- 0
  if(thinning <= 1) {
    post.bic[1] <- all.bic[1] <- bic.old
    post.model[1] <- model.old
    post.beta[1,] <- beta.old
    if (!is.null(B)) post.prob[1,] <- post.prob.old
    
    k <- 1
    if(verbose > 2)
      points(k, post.bic[k], cex = 0.5)
  }
  cont.accept <- numeric(0)
  beta.cont.accept <- rep(0,length(B))
  
  accept.fn <- function(x, move, fate = "accept") {
    move <- paste(fate, move, sep = ".")
    if(is.na(match(move, names(x))))
      x[move] <- 1
    else
      x[move] <- x[move] + 1
    x
  }
  
  
  ## Begin MCMC sampling
  for(i in 2:(nSamples*thinning+n.burnin)){
    ## Propose new network structure.
    if (rev.method=='node.edge'){
        M.new <-  propose.new.node.edge(M.old, max.parents,
        saved.scores = saved.scores, rev.method = rev.method,
        verbose = (verbose %in% c(2,4)), ...)  # Extension of Grzegorczyk and Husmeier (2008) Mach Learn 71: 265–305
        
        rev.ratio <- M.new$rev.ratio  ## When rev.method='nbhd', it has 1/2 prob of
                                      ## using Grzegorczyk and Husmeier (2008) method.
                                      ## In this case to reverse an edge, rev.ratio
                                      ## calculates the ratio between the proposal prob
                                      ## of reversing the edge, e.g. 1->2 and the proposal prob
                                      ## of reversing the edge, e.g. 2->1. Probability of choosing
                                      ## the edge is calculated later : ne.ratio <- sum(M.old)/sum(M.new)
                                      ## Otherwise, rev.ratio=1.
        
        move <- M.new$move ## Want to keep track of moves and acceptance rates.
        M.new <- M.new$M
        

        ## qtlnet:::nbhd.size() is different from Husmeier, D. (2003)'s neighbor size.
        ## For each node, it calculates possible downstream edge additons [node -> x], reversions [node -> x] or [x->node],
        ## and then, accumulates them.
        ## For deletions, it is the number of existing edges.
        ## Hence, qtlnet:::nbhd.size() counts reversions twice.
        ne.new <- qtlnet:::nbhd.size(M.new, max.parents)[[1]]
        
        ## Accept new model?
        ## Always if bic.new < bic.old.
        ## Otherwise do Metropolis-Hastings trick.
        if (move=='reverse-nbhd') ne.ratio <- sum(M.old)/sum(M.new) ## prob of choosing the edge
        else ne.ratio <- ne.old/ne.new ## single move's ratio

    }
    else if (rev.method %in% c('single','nbhd')){
        M.new <- propose.new.structure(M.old, max.parents,
        saved.scores = saved.scores, rev.method = rev.method,
        verbose = (verbose %in% c(2,4)), ...)  # If rev.method='single', Husmeier (2003) Bioinformatics 19:2271-2282. If rev.method='nbhd', Grzegorczyk and Husmeier (2008) Mach Learn 71: 265–305
        
        rev.ratio <- M.new$rev.ratio  ## When rev.method='nbhd', it has 1/2 prob of
                                      ## using Grzegorczyk and Husmeier (2008) method.
                                      ## In this case to reverse an edge, rev.ratio
                                      ## calculates the ratio between the proposal prob
                                      ## of reversing the edge, e.g. 1->2 and the proposal prob
                                      ## of reversing the edge, e.g. 2->1. Probability of choosing
                                      ## the edge is calculated later : ne.ratio <- sum(M.old)/sum(M.new)
                                      ## Otherwise, rev.ratio=1.
        
        move <- M.new$move ## Want to keep track of moves and acceptance rates.
        M.new <- M.new$M
        
        ## qtlnet:::nbhd.size() is different from Husmeier, D. (2003)'s neighbor size.
        ## For each node, it calculates possible downstream edge additons [node -> x], reversions [node -> x] or [x->node],
        ## and then, accumulates them.
        ## For deletions, it is the number of existing edges.
        ## Hence, qtlnet:::nbhd.size() counts reversions twice.
        ne.new <- qtlnet:::nbhd.size(M.new, max.parents)[[1]]

        
        ## Accept new model?
        ## Always if bic.new < bic.old.
        ## Otherwise do Metropolis-Hastings trick.
        if (move=='reverse-nbhd') ne.ratio <- sum(M.old)/sum(M.new) ## prob of choosing the edge
        else ne.ratio <- ne.old/ne.new ## single move's ratio
    }
    else if (rev.method %in% c('node','node.nbhd')){
        ## It picks a node and moves to the new structure by adding or reversing [x -> node], or by deleting [x -> node] or [node -> x].
        node <- sample(seq(n.pheno), 1)
        ## Gets possible changes.
        neighborM <- node.nbhd.size(M.old, node, max.parents)
        M.new <-  propose.new.node.structure(M.old, max.parents,
          saved.scores=saved.scores, node=node,  neighborM=neighborM, rev.method=rev.method, verbose=verbose)
       
        rev.ratio <- M.new$rev.ratio
        move <- M.new$move
        ne.new <- M.new$ne.new
        M.new <- M.new$M

        
        if (move=='reverse-nbhd') ne.ratio <- sum(M.old)/sum(M.new)
        else ne.ratio <- neighborM$nbhd.size/ne.new
    }
    
     
    ## Calculate P(G^{new}_Y|B, beta.old)
    # (prior probability of new structure given B and beta.old)
    if (!is.null(B)){
      Mprior.new <- score.prior(M.new,  beta.old, B)
    }
    
    ## Calculate the log likelihood of M.new.
    aux.new <- score.model.qtlnet.prior(M.new, saved.scores, cross, addcov, intcov,
                           threshold, nogenotype=nogenotype, verbose, method = method, ...)


    ## Typically only a few update.scores.
    tmp <- aux.new$update.scores
    if(!is.null(tmp)) {
      index <- nrow(saved.scores)
      index <- match(as.character(tmp$code), codes) + index * (tmp$pheno.col - 1)
      saved.scores[index] <- tmp$bic
    }

    bic.new <- aux.new$model.score
    model.new <- aux.new$model.name


    if(is.infinite(rev.ratio))
      mr <- 1
    else {
      mr <- exp(-0.5*(bic.new - bic.old)) * ne.ratio * rev.ratio
      if (!is.null(B)) mr <- mr * exp(Mprior.new - Mprior.old)
    }
    if(is.na(mr) | is.null(mr)){
        cat("mr:",mr,"\n")
        browser()
    }
    if(runif(1) <= min(1,mr)){
      M.old <- M.new
      bic.old <- bic.new
      ne.old <- ne.new
      model.old <- model.new
      if (!is.null(B)) post.prob.old <- bic.old - 2*Mprior.old + 2*partition.old
     if (!is.null(B)) Mprior.old <- Mprior.new
       
      cont.accept <- accept.fn(cont.accept, move)
    }
    else
      cont.accept <- accept.fn(cont.accept, move, "reject")

    
    if(!is.null(B)) {
      for (b in seq(B)){
        ## Sample new beta
        beta.new[b] <- propose.new.beta(beta.old[b], lambda[b])
        ## Calculate the partition function (normalizing term in the prior prob for phenotype networks)
        partition.new <- partition.eff(beta.new, n.pheno, B, count.parents.config, e.value)
        
        ## Calculate the prior probability of network in logarithm
        Mprior.new <- score.prior(M.old, beta.new, B)
        ## Calculate the prior probability of beta in logarithm
        prior.beta.new <- score.beta(beta.new, lambda)
        ## Calculate the ratio of prior probabilities
        mr <- exp(Mprior.new-Mprior.old + prior.beta.new-prior.beta.old + partition.old-partition.new)
        if(is.na(mr) | is.null(mr))
           browser()
       
       ## Accept the new beta with prob mr.
        if(runif(1) <= min(1, mr)){
          beta.old[b] <- beta.new[b]
          partition.old <- partition.new
          Mprior.old <- Mprior.new
          prior.beta.old <- prior.beta.new
          beta.cont.accept[b] <- beta.cont.accept[b] + 1
        }
      }  
    }

    

    ## Accumulate M's for average.
    if(i > n.burnin)
      Mav <- Mav + M.old
    
    ## Bookkeeping to save sample.
    if(((i-n.burnin)>0 & ((i-n.burnin) %% thinning) == 0)){      
      k <- k + 1
      all.bic[k] <- bic.new
      post.bic[k] <- bic.old
      post.model[k] <- model.old
      post.beta[k,] <- beta.old
      post.partition[k] <- partition.old
      if (!is.null(B)) post.prob[k] <- post.prob.old
     
      if(verbose) {
        print(c(i,k)) 
        if(verbose > 2)
          points(k, post.bic[k], cex = 0.5)
      }
      mav <- mav + M.old
      if (k %% 100 == 0) {
         out <- list(post.model = post.model,
              post.bic = post.bic, 
              Mav = Mav/(k*thinning), mav=mav/k,
              freq.accept = sum(cont.accept[grep('accept',names(cont.accept))]) / (k * thinning+n.burnin),
              cont.accept= cont.accept[sort(names(cont.accept))],
              saved.scores=saved.scores,
              all.bic=all.bic,
              cross = cross, post.beta=post.beta, post.partition=post.partition,
                  beta.freq.accept = beta.cont.accept/(k * thinning), post.prob=post.prob)  

            ## Attributes of qtlnet object.
          attr(out, "M0") <- M0
          attr(out, "threshold") <- threshold
          attr(out, "nSamples") <- nSamples
          attr(out, "intermediate.nSamples") <- k
          attr(out, "thinning") <- thinning
          attr(out, "pheno.col") <- pheno.col
          attr(out, "pheno.names") <- pheno.names
          attr(out, "addcov") <- addcov
          attr(out, "intcov") <- intcov
          attr(out, "burnin") <- burnin
          attr(out, "method") <- method
          attr(out, "random.seed") <- random.seed
          attr(out, "random.kind") <- RNGkind()
          attr(out, "B") <- B
          attr(out, "nogenotype") <- nogenotype
          attr(out, "lambda") <- lambda
  

         class(out) <- c("qtlnet","list","qtlnet.prior")
        
        save(out, file=filename, compress=TRUE)
      }
    }
  }
  

  ## Make average here.
  Mav <- Mav / (nSamples * thinning)
  cont.accept <- cont.accept[sort(names(cont.accept))]
  
  out <- list(post.model = post.model,
              post.bic = post.bic, 
              Mav = Mav, mav=mav/nSamples,
              freq.accept = sum(cont.accept[grep('accept',names(cont.accept))]) / (nSamples * thinning+n.burnin),
              cont.accept= cont.accept,
              saved.scores=saved.scores,
              all.bic=all.bic,
              cross = cross, post.beta=post.beta, post.partition=post.partition, beta.freq.accept = beta.cont.accept/(nSamples * thinning+n.burnin), post.prob=post.prob)

  ## Attributes of qtlnet object.
  attr(out, "M0") <- M0
  attr(out, "threshold") <- threshold
  attr(out, "nSamples") <- nSamples
  attr(out, "thinning") <- thinning
  attr(out, "pheno.col") <- pheno.col
  attr(out, "pheno.names") <- pheno.names
  attr(out, "addcov") <- addcov
  attr(out, "intcov") <- intcov
  attr(out, "burnin") <- burnin
  attr(out, "method") <- method
  attr(out, "random.seed") <- random.seed
  attr(out, "random.kind") <- RNGkind()
  attr(out, "B") <- B
  attr(out, "nogenotype") <- nogenotype
  attr(out, "lambda") <- lambda
  

  class(out) <- c("qtlnet","list","qtlnet.prior")
  
  out
}
################################################################
node.nbhd.size <- function(M, node, max.parents=3)
{
    
    ## This function finds a list of possible edges by deleting [x -> node] or [node->x],
    ## adding [x -> node], or reversing [x -> node]. 
    out <- NULL
    up.delete <- which(M[,node]==1)
    down.delete <- which(M[node,]==1)
    n.deletions <- length(c(up.delete, down.delete))
    if(length(up.delete)>0) {
        out <- data.frame(from=up.delete, to=node, move='delete')
        if (length(down.delete)>0){
            tmp <- data.frame(from=node, to=down.delete, move='delete')
            out <- rbind(out, tmp)
        }
    }
    else{
        if (length(down.delete)>0) out <- data.frame(from=node, to=down.delete, move='delete')
    }
    
    le <- ncol(M)
    up.additions <- upstream.additions(M, node, max.parents)
    n.additions <- length(up.additions)
    if (n.additions >0){
      if(!is.null(out)){
        tmp <- data.frame(from=up.additions, to=node, move='add')
        out <- rbind(out, tmp)
      } else {
        out <-data.frame(from=up.additions, to=node, move='add')
      }
    }
    
    n.reversions <- 0
    rev.allow <- check.upstream.reversions(M, node, max.parents)
    if (!is.null(rev.allow)) {
        n.reversions <- nrow(rev.allow)
        colnames(rev.allow) <- c('from','to')
        rev.allow <- data.frame(rev.allow)
        rev.allow <- cbind(rev.allow, move='reverse')
        if(!is.null(out)){
            out <- rbind(out, rev.allow)
        } else{
            out <-rev.allow
        }
    }
    
    nbhd.size <- n.deletions + n.additions + n.reversions
    out$move <- as.character(out$move)
    return(list(nbhd.size=nbhd.size, n.deletions=n.deletions,
    n.additions=n.additions, n.reversions=n.reversions, moves=out))
}

######################################################################
upstream.additions <- function(M, node, max.parents = 3)
{
    ## It finds a list of nodes that can be added (node -> x).
    
    
    ## Check for max.parents
    up <- which(M[,node]==1)
    if (length(up)>= max.parents) {
        ok.additions <- NULL
    }
    else{
      ## Forbidden upstream additions by making a cycle with downstreams
      upf <- which(M[node,] == 1)
      le <- length(upf)
      if (le>0) down <- check.downstream(M, upf)
      else down <- NULL
      forbidden <- unique(c(node, up, down))
      ok.additions <- seq(ncol(M))[-forbidden]
    }
    
    return(ok.additions)
}
################################################################
check.upstream.reversions <- function(M, node, max.parents = 3)
{
    ## Check on possible reversals: [node -> x] to [x -> node]
    ## allowed if no cycles produced.
    ## forbidden if cycles result.
    
    allowed <- NULL
    
    ## Check upstream.
    up <- which(M[,node] == 1)
    le <- length(up)
    
    if(le) {
        if(le == 1)
        allowed <- cbind(up, node)
        else { ## le > 1
            ## Multiple colliders.
            forbid.up <- rep(FALSE, le)
            for(k in 1:le)
            forbid.up[k] <- up[k] %in% check.upstream(M, up[-k])
            if(any(forbid.up)){
                if(any(!forbid.up))
                allowed <- cbind(up[!forbid.up], node)
            }
            else
            allowed <- cbind(up, node)
        }
    }
    
    ## Final check of max.parents.
    if(!is.null(allowed)) {
        wh <- which(apply(M[, allowed[,1], drop = FALSE], 2, sum) >= max.parents)
        if(length(wh)) {
            ## Remove pairs.
            allowed <- allowed[-wh,, drop = FALSE]
            if(nrow(allowed) == 0)
            allowed <- NULL
        }
    }
    
    allowed
}

#############################################
nbhd.size <- function(M, max.parents=3)
{
    ## Obsolte function.
    ## This function is slightly different from nbhd.size in qtlent.
    ## This function calculate possible number of structures from the current network structure
    ## by simple edge addition, deletion, and reversion
    n.deletions <- sum(M)
    n.additions <- 0
    n.reversions <- 0
    le <- ncol(M)
    wh.edge <- which(M==0)
    wh.diag<-diag(le)
    wh.diag <- which(wh.diag==1)
    wh.edge <- wh.edge[!(wh.edge %in% wh.diag)]
    for (w in wh.edge){
        n.additions <- n.additions + check.additions(M, w, max.parents)
    }
    
    wh.edge <- which(M==1)
    for(w in wh.edge){
        n.reversions <- n.reversions + check.reversions.edge(M,w,max.parents)
    }
    
    nbhd.size <- n.deletions + n.additions + n.reversions
    list(nbhd.size=nbhd.size,
    n.deletions=n.deletions,
    n.additions=n.additions,
    n.reversions=n.reversions)
}

######################################################
check.additions <- function(M, wh.edge, max.parents=3)
{
    ## check if adding the edge (wh.edge) in M
    ## still makes a DAG satisfying max.parents criteria
    ## Return: 1 for ok, 0 for no.
    M.new <- M
    M.new[wh.edge]<-1
    
    tmp <- (apply(M.new, 2, sum)>max.parents)
    if (sum(tmp)) return(0)
    
    check.cycle(M, wh.edge)
}

########################################################
check.reversions.edge <- function(M, wh.edge, max.parents=3)
{
    ## check if reversing the edge (wh.edge) in M
    ## still makes a DAG satisfying max.parents criteria
    ## Return: 1 for ok, 0 for no.
    le <- ncol(M)
    node.to <- ((wh.edge-1) %/% le) + 1
    node.from <- ((wh.edge-1) %% le) + 1
    
    M.new <- M
    M.new[node.from, node.to] <- 0
    M.new.tmp<-M.new
    M.new[node.to, node.from] <- 1
    
    tmp <- (apply(M.new, 2, sum)> max.parents)
    if (sum(tmp)) return(0)
    
    check.cycle(M.new.tmp, (node.from-1)*le+node.to)
}

#######################################################
check.cycle <- function(M, wh.edge){
    ## check if there is a cycle by adding wh.edge to M
    ## return value: 1 for without cycle, 0 for with cycles.
    le <- ncol(M)
    node.to <- ((wh.edge-1) %/% le) + 1
    node.from <- ((wh.edge-1) %% le) + 1
    
    is.down <- apply(M[node.to,,drop=FALSE], 2, sum)
    flag <- TRUE
    is.nocycle <- 1
    count.length <- 1
    while(flag){
        is.down<-which(is.down>0)
        if (length(is.down)==0) flag <- FALSE
        #flag <- length(is.down)
        if (any(is.down == node.from)){
            flag<-FALSE
            is.nocycle<-0
        }
        else {
            is.down <- apply(M[is.down,,drop=FALSE],2,sum)
            count.length <- count.length + 1
            if (count.length > le+1) flag <-FALSE
        }
    }
    return(is.nocycle)
}
################################################################
propose.new.beta <- function(beta.old, lambda)
{
  ## proposing a new beta
  ## proposal function : uniform distribution
  #runif(1, 0, 30)
  #b <- runif(1, beta.old - l, beta.old + l)
  #if (b > 10) b <- 20-b # 10-(b-10)
  #if (b>30) b<-60-b
  #if (b<0) b<- -1*b

  #b <- rexp(1, lambda)
  b <- runif(1, beta.old-1,beta.old+1)
  if (b<0) b<- -1*b
  b
  #rexp(1, 0.01)
}

##############################################################
score.prior <- function(M,  beta, B)
{
  ## Calculate the log of the numerator in the prior probability of a new model given beta
  if (!length(beta)) return(logical(0))
  tmp <- 0
  p<-dim(M)[1]
  for (i in seq(B)){
    #tmp<- tmp- beta.old[i]*sum(abs(M.new-B[[i]]))/(p*(p-1))
    #tmp<- tmp- beta.old[i]*(sum(abs(M.new-B[[i]]))-B.mean[[i]]*p)
    tmp<- tmp- beta[i]*sum(abs(M-B[[i]]))
  }
  tmp
}
#############################################################
score.beta <- function(beta, lambda)
{
  ## Caclualte the log of hyperprior for beta without a constant term
  ## P(beta | lambda) = lambda * exp(-lambda * beta)
  -1*sum(beta * lambda)
}

#############################################################
all.config <- function(n.pheno, max.parents.partition, B){
    ## Approximate calculation of node-wise expected value of |B-G| for each parent configuration and type of knowledge.
    
  if (max.parents.partition >= n.pheno) stop(paste("number of parents in partition function calculation by fan-in should be less than the number of phenotypes ", n.pheno))


  ## 1. number of possible parent combinations for each number of parents (1 to max.parents.partition)
  all.p <- sapply(seq(max.parents.partition), function(x) dim(combn(n.pheno-1, x))[2])
  ## total number of parent configurations
  all.p<-sum(all.p)
  all.parents.config <- vector("list", all.p)
  count<-1
  for(i in seq(max.parents.partition)){
    # each column is a combination
    x <- combn(n.pheno-1, i)
    for(xx in seq(dim(x)[2])){
      all.parents.config[[count]] <- x[,xx]
      count<-count+1
    }
  }


  ## 2. Expected value of |B-G| for each node, parent configuration, and type of knowledge.
  ## e.value is an array of (n.pheno, all.parents.config +1, number of types of knowledge).
  ## e.value[i,k,j] is a node-wise expected value approximation |B[[j]] - G_Y|
  ## for the child phenotype i, k-th parent configuration, and j-th knowledge,
  ## which is calculated by (1-B[[j]][k,i])+(B[[j]][-k,i]).
  
  ## +1 in "length(all.parents.config) + 1" for no parent case.
  e.value <- array(0, c(n.pheno, length(all.parents.config)+1, length(B)))
  for (i in seq(n.pheno)){
    for(b in seq(B)){
        e.value[i,1,b] <- sum(B[[b]][,i])  # no parent case for each phenotype
    }    
    for (j in seq(all.parents.config)){
      parents.config <- all.parents.config[[j]]
      parents.config[parents.config>=i] <- parents.config[parents.config>=i]+1
      for(b in seq(B)){
         #e.value[i,j,b]<- (sum(1-B[[b]][parents.config,i])+sum(B[[b]][-c(parents.config,i),i]) - B.mean[[b]])
        #e.value[i,j,b]<- sum(1-2*B[[b]][parents.config,i])
        e.value[i,j+1,b]<- (sum(1-B[[b]][parents.config,i])+sum(B[[b]][-c(parents.config,i),i]))
       }
    }
  }
  ## e.value: child-node-wise approximated expected value |B-G| for each parent configuration and type of knowledge
  ## count.parents.config: total number of parent configurations from no parent
  ##   up to max.parents.partition parents for a node.
  return(list(e.value=e.value, count.parents.config=length(all.parents.config)+1))
}


######################################################
partition.eff <- function(beta, n.pheno, B, count.parents.config, e.value){
    ## Calculate the log of partition function in an approximate way (normalizing term in
    ## prior probability of phenotype networks) specified by beta
    
    ## The parition function is defined to be \sum_{G_Y} exp(- \sum_j beta_j * |B[[j]] - G_Y|).
    ## Instead, it first calculates the node-wise |B-G| for each child node i approximately:
    ##  \sum_{k in node i's parent configuration} exp(- \sum_j beta_j * [(1-B[[j]][k,i]) + (B[[j]][-k,i])] ).
    ## Then, the product of them are used to approximate the partition function.
    ## it takes the product of the value for each node.
    
    ## 'e.value' is an array of (n.pheno, count.parents.config, number of types of knowledge).
    ## e.value[i,k,j] is a node-wise expected value approximation |B[[j]] - G_Y|
    ## for the child phenotype i, k-th parent configuration, and j-th knowledge,
    ## which is calculated by (1-B[[j]][k,i])+(B[[j]][-k,i]).
    
    ## 'count.parents.config' is the number of possible parents' configurations
    ## from no parent to max.parents.partition parents.
  Z <- 0
  for(i in seq(n.pheno)){
    node.exp <- 0
    for(j in seq(count.parents.config)){
      node.energy <- -1*sum(beta * e.value[i,j,])
      node.exp <- node.exp + exp(node.energy)
    }
    Z <- Z + log(node.exp)
  }
  #Z <- log(Z)
  Z
}



