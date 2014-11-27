##################################################################
## 2014/09/16
## These functions are obtained from an R package : qtlnet (version 1.2.5
## Copyright (C) 2012 Elias Chaibub Neto and Brian S. Yandell) from
## CRAN repository. They are redistributed under the terms of the GNU 
## General Public License as published by the Free Software Foundation.
##
## Note: propose.new.structure is modified a little to use a mix of
## 'single' method and 'nbhd' method when rev.method='nbhd' to reverse an edge
##
## Functions: propose.add, propose.new.node.edge, propose.new.structure,
##     update.node, drop.edge, sample.exclude.down, reverse.proposal,
##     rev.edge, find.index.parent, index.parents, find.parent.score, 
##     forbidden.additions, check.qtlnet, check.upstream, check.downstream,
##     check.reversions, print.qtlnet
## New function: propose.new.node.structure
######################################################################
propose.add <- function(M, node, max.parents)
{
  ## This function samples a child node from the given node.
    
  ## Add an edge with a direction.
  aux1 <- forbidden.additions(M, node, max.parents)
  ## Is there a down option?
  aux3 <- unique(c(node, aux1$upf, aux1$downf))
  aux3 <- seq(ncol(M))[-aux3]
  
  if(length(aux3)) {
    ## Check if any of aux3 is at or above max.parents.

    wh <- which(apply(M[, aux3, drop = FALSE], 2, sum) >= max.parents)
    
    if(length(wh))
      aux3 <- aux3[-wh]
  }
  if(length(aux3) == 0)
    aux3 <- 0
  else {
    if(length(aux3) > 1)
      aux3 <- sample(aux3, 1)
  }
  aux3
}  
######################################################################
propose.new.node.edge <- function(M, max.parents = 3,
                                  saved.scores, rev.method = "node.edge",
                                  propose = c(node = 1, edge = 2, reverse = 10, drop = 2),
                                  check.down = FALSE,
                                  verbose = FALSE, ...)
{
  ## Ref: Grzegorczyk and Husmeier (2008) Mach Learn 71: 265–305.
  ## Extension is to updating individual nodes and dropping edges.
  rev.ratio <- 1

  ## Proposal weights for types of changes.
  if(any(propose <= 0))
    stop("propose values must be positive")
  name.propose <- names(propose)
  propose <- array(propose, 4)
  if(length(name.propose) == 4)
    names(propose) <- name.propose
  else
    names(propose) <- c("node","edge","reverse","drop")
  
  ## To speed up, use M as logical matrix. Still return 0/1 matrix for now.
  le.nodes <- ncol(M)
  le.edges <- sum(M)

  prob.node <- propose["node"] * le.nodes
  prob.node <- prob.node / (prob.node + propose["edge"] * le.edges)
  if(runif(1) <= prob.node) {
    move <- "node"
    ## Propose new parents for a node.
    node <- sample(seq(le.nodes), 1)
    if(verbose)
      cat(move, "")
    
    aux2 <- update.node(M, max.parents, saved.scores, node, verbose)
    rev.ratio <- aux2$rev.ratio
    M <- aux2$M
    if(verbose) {
      up <- qtlnet:::node.parents(M, node)$parents
      cat(node, "up:", up, "\n")
    }
  }
  else {
    ## Propose change to edge.
    move <- "edge"
    if(verbose)
      cat(move, "")

    ## Pick a valid edge.
    edge <- which(M == 1)
    if(length(edge) > 1)
      edge <- sample(edge, 1)
    edge <- c(row(M)[edge], col(M)[edge])
    if(verbose)
      cat("edge ")
      
    ## Cyclic mistake checks.
    ## These should not happen!
    if(edge[1] == edge[2]) {
      cat("edge identity:", edge, "\n")
      browser()
    }
    aux1 <- check.downstream(M, edge[2])
    if(any(aux1 == edge[1])) {
      cat("edge downfall:", edge, aux1, "\n")
      browser()
    }
    
    prob.reverse <- propose["reverse"]
    prob.reverse <- prob.reverse / (prob.reverse + propose["drop"])
    if(runif(1) <= prob.reverse) {
      ## Propose to reverse an edge.
      if(verbose)  cat("reverse ")
      aux2 <- rev.edge(M, max.parents, saved.scores, edge, verbose)
    }
    else {
      ## Propose to drop an edge.
      if(verbose)  cat("drop ")
      aux2 <- drop.edge(M, max.parents, saved.scores, edge, verbose)
    }
    rev.ratio <- aux2$rev.ratio
    M <- aux2$M
    
    if(verbose) {
      for(i in 1:2) {
        up <- qtlnet:::node.parents(M, edge[i])$parents
        cat(edge[i], "up:", up, "\n")
      }
    }

    if(check.down) {
      ## NOTE: This takes time and should be dropped eventually.
      ## Check if we somehow have a created a cycle earlier. This should not happen.
      down <- check.downstream(M, edge[2])[-1]
      up <- check.upstream(M, edge[2])[-1]
      if(length(down) & length(up)) {
        if(any(down %in% up)) {
          cat("downup:", edge, "\n", sort(down), "\n", sort(up), "\n")
          browser()
        }
      }
    }
  }

  
  list(M = M, rev.ratio = rev.ratio, move = move)
}
######################################################################
propose.new.node.structure <- function(M, max.parents = 3,
saved.scores, node, neighborM, rev.method = c("node","node.nbhd"),
verbose = FALSE)
{
    
    ## It samples a node and proposes a move by adding or reversing [x -> node],
    ## or by deleting [x -> node] or [node -> x].
    ## This code is modified from qtlnet:::propose.new.structure().

    rev.ratio <- 1
    ne.new <- 0
    ## To speed up, use M as logical matrix. Still return 0/1 matrix for now.
    le.nodes <- ncol(M)

    flag <- TRUE
    
    if (rev.method=='node.nbhd') move <- sample(c('reverse-nbhd','single'), 1)
    else move<-'single'
    
    if (move=='reverse-nbhd'){
        ## Reverse method of Grzegorczyk and Husmeier.
        aux1 <- c(which(M[,node] == 1), le.nodes + which(M[node,] == 1))
        if(!(flag <- !length(aux1))) {
            if(length(aux1) > 1) aux1 <- sample(aux1, 1)
            if(aux1 > le.nodes) aux3 <- c(node, aux1 - le.nodes)
            else aux3 <- c(aux1, node)
            
            aux1 <- check.downstream(M, aux3[2])
            if(any(aux1 == aux3[1])) {
                ## This should not happen!
                cat(move, "downfall:", aux3[1], aux1, "\n")
                browser()
            }
            
        }
        if(!flag) {
            if(aux3[1] == aux3[2]) {
                ## This should not happen!
                cat(move, "identity:", aux3, "\n")
                browser()
            }
            aux2 <- rev.edge(M, max.parents, saved.scores, aux3, verbose)
            rev.ratio <- aux2$rev.ratio
            if (is.na(rev.ratio)){
                cat('rev.ratio is NA in reverse-nbhd move\n')
                browser()
            }
            M <- aux2$M
        }
    }
    
    ## Single move by changing [x -> node] by add, delete, or reverse from node.nbhd.size()
    if (flag){
       aux1 <- sample(seq(neighborM$nbhd.size), 1)
       move <- neighborM$moves[aux1,'move']
       aux3 <- neighborM$moves[aux1,c(1,2)]
       switch(move,
       add = {
          M[aux3[1,1],aux3[1,2]] <- 1
       },
       reverse = {
         M[aux3[1,1],aux3[1,2]] <- 0
         M[aux3[1,2],aux3[1,1]] <- 1
       },
       delete = {
         if (aux3[1,1]==node) M[aux3[1,1],aux3[1,2]] <- 0
         else M[aux3[1,2],aux3[1,1]] <- 0
       })
       ## proposal probability from the new one to the old one.
       if (move=='add'){
           ## there are two nodes possibly selected to from M.new to M.old
           ne.new <- node.nbhd.size(M, aux3[1,1], max.parents)$nbhd.size
           ne.new <- ne.new/2 + neighborM$nbhd.size/2
       }
       else ne.new <- node.nbhd.size(M, aux3[1,1], max.parents)$nbhd.size
       if (ne.new==0) {
           cat('no move possible.\n')
           browser()
       }
    }
    if (verbose) cat(move,"\n")
    list(M = M, rev.ratio = rev.ratio, move = move, ne.new = ne.new)
}

######################################################################
propose.new.structure <- function(M, max.parents = 3,
                                  saved.scores, rev.method = "nbhd",
                                  verbose = FALSE, ...)
{
  ## rev.method ='single': Husmeier, D. (2003) Bioinformatics 19:2271-2282
  ## Acceptance rate is <20%. Could we improve here?
  ## Yes. rev.method='nbhd': it is an update version of 'rev.method=single'
  ## method by using reversing the edge method described in
  ## Grzegorczyk and Husmeier (2008) Mach Learn 71: 265–305.
  
  ## JY thinks that rev.method='single' is not exactly what Husmeier (2003) does.
  ## It is because in the code below, when it selects a node and its move,
  ## the prob is not the same. For example, 1->2->3,
  ## Node selected   |  move                | new structure | prob
  ##   1             |  add 1->3            | (1)(2|1)(3|1,2)| 1/3 *1
  ##   1             |  delete 1->2         | (1)(2)(3|2)   | 1/3 * 1
  ##   1             |  reverse 1->2        | (1|2)(2)(3|2) | 1/3 * 1
  ##   2             |  delete 1->2         | (1)(2)(3|2)   | 1/3*1/2
  ##   2             |  delete 2->3         | (1)(2|1)(3)   | 1/3*1/2
  ##   2             |  reverse 1->2        | (1|2)(2)(3|2) | 1/3*1/2
  ##   2             |  reverse 2->3        | (1)(2|1,3)(3) | 1/3*1/2
  ##   3             |  delete 2->3         | (1)(2|1)(3)   | 1/3*1
  ##   3             |  reverse 2->3        | (1)(2|1,3)(3) | 1/3*1
  
  ## new structure  | total proportion
  ## (1)(2|1)(3|1,2) | 2/6
  ## (1)(2)(3|2)     | 3/6
  ## (1|2)(2)(3|2)   | 3/6
  ## (1)(2|1)(3)     | 3/6
  ## (1)(2|1,3)(3)   | 3/6
  
  ## qtlnet:::nbhd.size() is different from Husmeier, D. (2003)'s neighbor size.
  ## Also, it is different from what this function's proposal probability.
  
  rev.ratio <- 1

  ## To speed up, use M as logical matrix. Still return 0/1 matrix for now.
  le.nodes <- ncol(M)
  flag <- TRUE
  while(flag){
    ## Pick node and decide on add/delete/reverse.
    ## Keep doing this until successful.
    
    node <- sample(seq(le.nodes),1)
    #move <- sample(c("add","delete","reverse"),1, prob=c(0.2,0.4,0.6)) JY changed.
    move <- sample(c("add","delete","reverse"),1, prob=c(0.2,0.2,0.6))
    if(verbose)
      cat(node, move, "")

    switch(move,
           add = {
               aux3 <- propose.add(M, node, max.parents)  # samples a child node from the given node
             if(!(flag <- (aux3 == 0))) {
               aux1 <- check.downstream(M, aux3)
               if(any(aux1 == node)) {
                 ## This should not happen!
                 cat(move, "downfall:", node, aux1, "\n")
                 browser()
               }
               M[node,aux3] <- 1
               if(verbose)
                 cat(aux3, "\n")
             } 
           },
           reverse = {
             ##JY modifying. MCMC uses Grzegorczyk and Husmeier (2008)'s method with prob 1/2
             ## and single method with prob 1/2
             rev.method2 = sample(c(rev.method,'single'), size=1, prob=c(0.2,0.4))
             
             ## Reverse direction of an edge.
             #if(flag <- (rev.method == "single")) {
             if(flag <- (rev.method2 == "single")) {            
               aux1 <- check.reversions(M, node, max.parents)
               if(!(flag <- is.null(aux1))) {
                 le.rev <- nrow(aux1)
                 aux3 <- sample(seq(le.rev), 1)
                 aux3 <- aux1[aux3, ]
                 M[aux3[1], aux3[2]] <- 0
                 M[aux3[2], aux3[1]] <- 1
               }
             }
             else { ## Reverse method of Grzegorczyk and Husmeier.
               aux1 <- c(which(M[,node] == 1), le.nodes + which(M[node,] == 1))
               if(!(flag <- !length(aux1))) {
                 if(length(aux1) > 1)
                   aux1 <- sample(aux1, 1)
                 if(aux1 > le.nodes)
                   aux3 <- c(node, aux1 - le.nodes)
                 else
                   aux3 <- c(aux1, node)

                 aux1 <- check.downstream(M, aux3[2])
                 if(any(aux1 == aux3[1])) {
                   ## This should not happen!
                   cat(move, "downfall:", aux3[1], aux1, "\n")
                   browser()
                 }
                 
               }
               if(!flag) {
                 if(aux3[1] == aux3[2]) {
                   ## This should not happen!
                   cat(move, "identity:", aux3, "\n")
                   browser()
                 }
                 aux2 <- rev.edge(M, max.parents, saved.scores, aux3, verbose)
                 rev.ratio <- aux2$rev.ratio
                 M <- aux2$M
               }
             }
             if(!flag & verbose)
                 cat(aux3[1], aux3[2], "\n")
           },
           delete = {
             ## Delete an existing edge through node.
             aux1 <- c(which(M[, node] == 1), le.nodes + which(M[node,] == 1))
             if(!(flag <- !length(aux1))) {
               if(length(aux1) > 1)
                 aux1 <- sample(aux1, 1)
               if(aux1 > le.nodes)
                 aux3 <- c(node, aux1 - le.nodes)
               else
                 aux3 <- c(aux1, node)
               M[aux3[1], aux3[2]] <- 0
               if(verbose)
                 cat(aux3, "\n")
             }
           })
  }

  ## Check if we somehow have a created a cycle earlier. This should not happen.
  down <- check.downstream(M, aux3[2])[-1]
  up <- check.upstream(M, aux3[2])[-1]
  if(length(down) & length(up)) {
    if(any(down %in% up)) {
      cat(move, "downup:", aux3, "\n", sort(down), "\n", sort(up), "\n")
      browser()
    }
  }

  if (move=='reverse') move <- paste(move, rev.method2, sep='-')
  list(M = M, rev.ratio = rev.ratio, move = move)
}
######################################################################
update.node <- function(M, max.parents, saved.scores, node, verbose = FALSE)
{
  ## Scores for possible parent sets.
  down <- check.downstream(M, node)[-1]
  index <- index.parents(saved.scores, down, node)
  if(length(index))
    z.1 <- saved.scores[-index, node]
  else 
    z.1 <- saved.scores[, node]
  z.1 <- exp(min(z.1) - z.1)
  s.1 <- sum(z.1)
  if(length(z.1) == 0) ## Should not happen
    browser()
    
  ## Find probability for current parents of node.
  rev.ratio <- z.1[find.parent.score(M, saved.scores, -index, node)]
    
  ## Make node an orphan.
  M[,node] <- 0

  ## Sample new parents of node.
  parent <- sample(seq(length(z.1)), 1, prob = z.1)

  ## Proposal ratio.
  rev.ratio <- z.1[parent] / rev.ratio

  ## Add edges to graph for new parents of node.
  new.parent <- find.index.parent(M, saved.scores, -index, node, parent)
  if(length(new.parent))
    M[new.parent, node] <- 1

  if(verbose) {
    cat("node rev.ratio", rev.ratio, "\n")
    if(is.na(rev.ratio))
      browser()
  }
  
  list(M = M, rev.ratio = rev.ratio)

}
######################################################################
drop.edge <- function(M, max.parents, saved.scores, node.pair, verbose = FALSE)
{
  ## Grzegorczyk and Husmeier (2008) Mach Learn 71: 265–305.
  ## Call provides node.pair[1 -> 2].

  ## Step 1: Orphan both nodes after computing reverse proposal prob.
  rev.ratio <- reverse.proposal(M, saved.scores, node.pair, verbose)
  M[, node.pair] <- 0

  ## Step 2. Sample new parents for node 1 without node 2 as one parent.
  
  ## Exclude parents downstream of 1 or 2.
  tmp <- sample.exclude.down(M, saved.scores, node.pair, verbose)
  q.1 <- tmp$prob
  parent <- tmp$parent
  
  ## Add edges to graph for parents of node 1.
  if(length(parent))
    M[parent, node.pair[1]] <- 1

  ## Step 3. Sample new parents for node 2 without node 1 as parent.

  ## Exclude parents downstream of 1 or 2.
  tmp <- sample.exclude.down(M, saved.scores, rev(node.pair), verbose)
  q.2 <- tmp$prob
  parent <- tmp$parent

  ## Proposal ratio.
  rev.ratio <- q.1 * q.2 / rev.ratio

  ## Add edges to graph for parents of node 2.
  if(length(parent))
    M[parent, node.pair[2]] <- 1

  if(verbose) {
    cat("edge rev.ratio", rev.ratio, "\n")
    if(is.na(rev.ratio))
      browser()
  }
  
  list(M = M, rev.ratio = rev.ratio)
}
######################################################################
sample.exclude.down <- function(M, saved.scores, node.pair, verbose = FALSE)
{
  ## Exclude parents downstream of 1 or 2.
  down <- unique(c(check.downstream(M, node.pair[2]),
                   check.downstream(M, node.pair[1])[-1]))
  index <- index.parents(saved.scores, down, node.pair[1])
  z.1 <- saved.scores[-index, node.pair[1]]
  z.1 <- exp(min(z.1) - z.1)
  if(length(z.1) == 0) ## Should not happen
    browser()
  parent <- sample(seq(length(z.1)), 1, prob = z.1)
  q.1 <- z.1[parent] / sum(z.1)

  new.parent <- find.index.parent(M, saved.scores, -index, node.pair[1], parent)

  list(parent = new.parent, prob = q.1)
}
######################################################################
reverse.proposal <- function(M, saved.scores, node.pair, verbose = FALSE)
{
  ## Parents of 2 include 1 but exclude all downstream of 2.
  index <- index.parents(saved.scores, node.pair[1], node.pair[2])
  down <- check.downstream(M, node.pair[2])[-1]
  if(length(down))
    index <- index[-index.parents(saved.scores[index,, drop = FALSE], down, node.pair[2])]
  z.2 <- saved.scores[index, node.pair[2]]
  z.2 <- exp(min(z.2) - z.2)
  q.2 <- z.2[find.parent.score(M, saved.scores, index, node.pair[2])] / sum(z.2)

  ## Parents of 1 exclude nodes at or downstream of 2.
  down <- unique(c(check.downstream(M, node.pair[2]),
                   check.downstream(M, node.pair[1])[-1]))
  index <- index.parents(saved.scores, down, node.pair[1])
  z.1 <- saved.scores[-index, node.pair[1]]
  z.1 <- exp(min(z.1) - z.1)
  q.1 <- z.1[find.parent.score(M, saved.scores, -index, node.pair[1])] / sum(z.1)

  rev.ratio <- q.1 * q.2

  if(verbose) {
    cat("rev.ratio", rev.ratio, "\n")
    if(is.na(rev.ratio))
      browser()
  }
  rev.ratio
}

######################################################################
rev.edge <- function(M, max.parents, saved.scores, node.pair, verbose = FALSE)
{
  ## Grzegorczyk and Husmeier (2008) Mach Learn 71: 265–305.
  ## Call provides node.pair[1 -> 2].

  ## Step 1: Orphan both nodes after computing reverse proposal prob.
  rev.ratio <- reverse.proposal(M, saved.scores, node.pair, verbose)
  M[, node.pair] <- 0

  ## Step 2. Sample new parents for node 1 with node 2 as one parent.
  
  ## Exclude parents downstream of 1 (except 2).
  index <- index.parents(saved.scores, node.pair[2], node.pair[1])
  down <- check.downstream(M, node.pair[1])[-1]
  if(length(down))
    index <- index[-index.parents(saved.scores[index,, drop = FALSE], down, node.pair[1])]
  z.1 <- saved.scores[index, node.pair[1]]
  z.1 <- exp(min(z.1) - z.1)
  if(length(z.1) == 0) ## Should not happen
    browser()
  parent <- sample(seq(length(z.1)), 1, prob = z.1)
  q.1 <- z.1[parent] / sum(z.1)
  new.parent <- find.index.parent(M, saved.scores, index, node.pair[1], parent)

  ## Add edges to graph for parents of node 1.
  if(length(new.parent))
    M[new.parent, node.pair[1]] <- 1

  ## Step 3. Sample new parents for node 2 without node 1 as parent.

  ## Exclude parents downstream of 1 or 2.
  ## Somehow this is leading to cycles!
  ## SOmehow the rev(node.pair) is not working.
  tmp <- sample.exclude.down(M, saved.scores, rev(node.pair))
  q.2 <- tmp$prob
  new.parent <- tmp$parent

  ## Add edges to graph for parents of node 2.
  if(length(new.parent))
    M[new.parent, node.pair[2]] <- 1

  rev.ratio <- q.1 * q.2 / rev.ratio

  if(verbose) {
    cat("edge rev.ratio", rev.ratio, "\n")
    if(is.na(rev.ratio))
      browser()
  }
  
  list(M = M, rev.ratio = rev.ratio)
}
######################################################################
find.index.parent <- function(M, saved.scores, index, node, new.parent)
{
  if(length(index))
    parents <- dimnames(saved.scores)[[1]][index][new.parent]
  else
    parents <- dimnames(saved.scores)[[1]][new.parent]
  unlist(sapply(strsplit(parents, ",", fixed = TRUE),
                function(x, node) {
                  x <- as.numeric(x)
                  x + (x >= node)
                }, node))
}
######################################################################
index.parents <- function(saved.scores, node1, node2)
{
  ## Row names of saved.scores renumber around missing node.
  parents <- node1 - (node2 < node1)

  ## Which parents include node.pair[2]?
  which(sapply(strsplit(dimnames(saved.scores)[[1]], ",", fixed = TRUE),
               function(x, parents) any(x %in% parents),
               parents))
}
######################################################################
find.parent.score <- function(M, saved.scores, index, node)
{
  parents <- qtlnet:::node.parents(M, node)$parents
  if(!is.null(parents)) {
    parents <- parents - (parents > node)
    if(length(index))
      match(paste(parents, collapse = ","), dimnames(saved.scores)[[1]][index])
    else
      match(paste(parents, collapse = ","), dimnames(saved.scores)[[1]])
  }
  else
    1
}
######################################################################
forbidden.additions <- function(M, node, max.parents = 3)
{
  ## It finds a list of nodes that cannot be added as a child node from the given node.
  ## The list consists of current child nodes for the node in M (upf) and
  ## current parents or ancestors for the node in M (downf).
  
  ## upf are forbidden upstream nodes (already present downstream).
  ## downf are forbidden downstream nodes (already present upstream).
  ## Check on cycles is one step. See check.reversions() for longer cycles.

  ## Forbidden upstream additions
  ## Forbidden upstream additions are direct child of the node (existing edge from the node)
  upf <- which(M[node,] == 1)
  if(length(upf)==0)
    upf <- NULL

  ## Forbidden downstream additions. More complicated.
  ## Forbidden downstreams additions are checked for cycles.
  downf <- which(M[,node] == 1)
  le <- length(downf)
  
  ## If at least max.parents are causal for node, forbid any more.
  if(le >= max.parents) 
    downf <- seq(nrow(M))[-node]
  else {
    if(le > 0)
      downf <- check.upstream(M, downf)
    else
      downf <- NULL
  }
  
  list(upf=upf, downf=downf)
} 
######################################################################
check.qtlnet <- function(object,
                         min.prob = 0.9,
                         correct = TRUE,
                         verbose = FALSE,
                         ...)
{
  pheno.names <- attr(object, "pheno.names")
  n.pheno <- length(pheno.names)

  forbid <- 1
  while(!is.null(forbid)) {
    M <- threshold.net(object, min.prob = min.prob, ...)
    min.prob <- attr(M, "min.prob")
    M1 <- 1 * (M > 0)

    ## This may not be right yet.
    forbid <- NULL
    for(i in seq(n.pheno)) {
      downf <- check.upstream(M1, i)
      wh <- which(M1[i, downf] == 1)
      if(length(wh))
        forbid <- rbind(forbid, cbind(downf[wh],i, M[i, downf[wh]]))
    }
    if(!is.null(forbid)) {
      forbid <- data.frame(forbid)
      names(forbid) <- c("cause","react","prob")
      for(i in 1:2)
        forbid[[i]] <- ordered(pheno.names[forbid[[i]]], pheno.names)
    }
    if(verbose)
      print(forbid)
  }
  attr(M, "min.prob") <- min.prob
  
  list(forbid = forbid, M = M)
}
######################################################################
check.upstream <- function(M, nodes)
{
  ## After accounting for scanone, 85% of time is nbhd.size.
  ## Of that, 85% (75% overall) is in check.upstream.

  count.upok <- nrow(M) - length(nodes)
  
  ## check upstream to see if a cycle would be created.
  is.up <- apply(M[, nodes, drop = FALSE], 1, sum)
  flag <- TRUE
  while(flag){
    aux1 <- which(is.up > 0)
    if(flag <- (length(aux1) > 0 & count.upok > 0)) {
      aux1 <- aux1[!(aux1 %in% nodes)]
      if(flag <- length(aux1)) {
        is.up <- apply(M[, aux1, drop = FALSE], 1, sum)
        nodes <- c(nodes, aux1)
        count.upok <- count.upok - flag
      }
    }
  }
  nodes
}
######################################################################
check.downstream <- function(M, nodes)
{
  count.downok <- nrow(M) - length(nodes)
  
  ## check downstream to see if a cycle would be created.
  is.down <- apply(M[nodes,, drop = FALSE], 2, sum)
  flag <- TRUE
  while(flag){
    aux1 <- which(is.down > 0)
    if(flag <- (length(aux1) > 0 & count.downok > 0)) {
      aux1 <- aux1[!(aux1 %in% nodes)]
      if(flag <- length(aux1)) {
        is.down <- apply(M[aux1,, drop = FALSE], 2, sum)
        nodes <- c(nodes, aux1)
        count.downok <- count.downok - flag
      }
    }
  }
  nodes
}
######################################################################
check.reversions <- function(M, node, max.parents = 3)
{
  ## Check on possible reversals of directions.
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
    
  ## Check downstream.
  down <- which(M[node,] == 1)
  le <- length(down)
  
  allowed2 <- NULL
  
  if(le) {
    if(le == 1)
      allowed2 <- cbind(node, down)
    else { ## le > 1
      ## Multiple colliders.
      forbid.down <- rep(FALSE, le)
      for(k in 1:le)
        forbid.down[k] <- down[k] %in% check.downstream(M, down[-k])
      if(any(forbid.down)){
        if(any(!forbid.down))
          allowed2 <- cbind(node, down[!forbid.down])
      }
      else
        allowed2 <- cbind(down, node)
    }
  }
  allowed <- rbind(allowed, allowed2)

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
print.qtlnet <- function (x, cutoff = 0.01, digits = 3, ...)
{
    cat("\nModel averaged probabilities for edge direction (row -> col):\n")
    mav <- qtlnet:::get.model.average(x)
    print(round(mav, digits))
    cat("\nPosterior probabilities by causal model:\n")
    pp <- qtlnet:::get.posterior.prob(x)
    cutoff <- min(cutoff, max(pp$post.prob))
    wh <- which(pp$post.prob >= cutoff)
    print(pp[wh, , drop = FALSE])
    invisible(list(mav = mav, pp = pp))
}

