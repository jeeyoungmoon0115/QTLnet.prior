#############################################################
## 2014/07/17
##  Copyright (C) 2014 Jee Young Moon and Brian S. Yandell
## This code is under the GNU.
##
## Note:
##   score.model.qtlnet.prior() is modified qtlnet:::score.model()
##      from in R/qtlnet, to take care of no QTL mapping case.
##   calc.pheno.bic() is modified from qtlnet:::calc.bic() in R/qtlnet, which does no QTL mapping.
##   scan.pheno() is modified from qtlnet:::scan.genome() in R/qtlnet, which does no QTL mapping.
##   bic.qtlnet.pheno() is modified from qtlnet:bic.qtlnet() in R/qtlnet, which does no QTL mapping.
##
## Functions: score.model.qtlnet.prior, calc.pheno.bic, scan.pheno,
##        bic.qtlnet.pheno
##############################################################

score.model.qtlnet.prior <- function(M, saved.scores, cross, addcov, intcov, threshold, nogenotype,
                        verbose = TRUE, ...)
{
    ## This function calculates the log-likelihood of network M.
    ## This function is modified from qtlnet:::score.model() to consider no QTL mapping case.
  n.pheno <- ncol(M)
  model.score <- 0
  mymodel <- rep(NA,n.pheno)
  count.score <- 0
  update.scores <- list(code = NULL, pheno.col = NULL, bic = NULL)
  
  for(i in 1:n.pheno){
    pheno <- qtlnet:::node.parents(M, i)
    mymodel[i] <- paste("(",paste(pheno$identifier,")",sep=""),sep="")

    ## Find saved BIC score if already computed.
    bic <- qtlnet:::find.bic(pheno$code, i, update.scores, saved.scores, ...)

    if(is.na(bic)){
      if (nogenotype){
        #### QTL mapping is not done. calc.pheno.bic instead of qtlnet:::calc.bic
        run <- calc.pheno.bic(cross, pheno$code, i, pheno$parents,
                      addcov, intcov, threshold,
                      n.pheno, update.scores = update.scores,
                      saved.scores = saved.scores, ...)
      }
      
      if (!nogenotype){
          run <- qtlnet:::calc.bic(cross, pheno$code, i, pheno$parents,
          addcov, intcov, threshold, n.pheno, update.scores = update.scores,
          saved.scores = saved.scores, ...)
      }
      
      ## Update saved scores.
      update.scores$code <- c(update.scores$code, as.character(run$code))
      update.scores$pheno.col <- c(update.scores$pheno.col, run$pheno.col)
      update.scores$bic <- c(update.scores$bic, run$bic)

      ## Print scan info if verbose.
      count.score <- count.score + 1
      if(verbose) {
        if(count.score == 1) cat("\nscan ")
        cat(paste("(", paste(run[, "pheno.col"], collapse = ","), "|",
                  paste(pheno$parents, collapse = ","), ")", sep = ""), "")
      }
      bic <- run[i == run[, "pheno.col"], "bic"]
    }
    ## Accumulate model score.
    model.score <- model.score + bic
  }
    
  list(model.score = model.score,
       update.scores = update.scores,
       model.name = paste(mymodel,collapse=""))
}
###########################################################################################
#agree.covs <- function(x,y) {
#  out <- length(x) == length(y)
#  if(out & length(x))
#    out <- all(sort(x) == sort(y))
#  out
#}
###########################################################################################
#find.bic <- function(code, pheno.col, update.scores = NULL, saved.scores #= NULL, ...)
#{
#  bic <- rep(NA, length(code))
#  
#  if(!is.null(saved.scores)) {
#    n.scores <- nrow(saved.scores)
#    score.pointer <- match(code, dimnames(saved.scores)[[1]])
#    bic <- saved.scores[score.pointer + n.scores * (pheno.col - 1)]
#  }
#  bic.na <- is.na(bic)
#  if(any(bic.na) & !is.null(update.scores)) {
#    wh <- match(paste(code[bic.na], pheno.col[bic.na], sep = "."),
#                paste(update.scores$code, update.scores$pheno.col, sep = #"."),
#                nomatch = 0)
#    if(!all(wh == 0))
#      bic[bic.na][wh > 0] <- update.scores$bic[wh]
#  }
#  bic 
#}
######################################################################
#node.parents <- function(M, node)
#{
#  aux <- which(M[,node] == 1)
#  if(length(aux) == 0){ 
#    parents <- NULL
#    identifier <- as.character(node)
#  }
# else{
#    parents <- aux
#    identifier <- paste(node,paste(parents,collapse=","),sep="|")
#  }
#  code <- find.code(node, parents)
#  list(parents=parents, identifier=identifier, code = code)
#}
#####################################################################
calc.pheno.bic <- function(cross, code, pheno.col, parents, addcov, intcov, threshold,
                     n.pheno, method = "hk", min.nonmissing = 50, ...)
{
    ## It does not run QTL mapping. This is for the case when there
    ## is no genotype data.
  
  ## Match phenotypes with same parents to run scan across genome.
  run <- qtlnet:::match.parents(cross, code, pheno.col, parents, addcov, intcov,
                       n.pheno, ...)
  
  ## Run pheno.cols through scanone in scan.genome.
  if(sum(!run$na) < min.nonmissing) {
    warning(paste("BIC set to 0: only ", sum(!run$na),
                  " non-missing values for combination: (",
                  paste(run$pheno.cols, collapse = ","), "|",
                  paste(parents, collapse = ","), ")", sep = ""))
    bic <- rep(0, length(run$pheno.cols))
  }
  else ## QTL mapping is not done. scan.pheno instead of scan.genome
    bic <- scan.pheno(subset(cross, ind = !run$na),
                       run$pheno.cols, parents,
                       addcov[[pheno.col]], intcov[[pheno.col]],
                       threshold[run$pheno.cols], method, ...)
  run$bic <- bic
  data.frame(code = run$code,
             pheno.col = run$pheno.cols,
             bic = bic)
}
#####################################################################
scan.pheno <- function(cross, pheno.col, pheno.parents, addcov, intcov,
                        threshold, method, ...)
{
  ## Modified version of scan.genome
  ## It does not QTL mapping. This is for the case when there
  ## is no genotype data.

  ## This currently scans one phenotype at a time.
  ## Want to extend to all phenotypes not in names(covM.dat).
  ## To do that, assume covariates are same for all i.
  
  n.pheno <- length(pheno.col)

  ## *** The following could be simplified.
  ## *** However, it involves a rethinking of myformula.
  ## *** Low priority for now.
  
  ## Design matrices for qtl/scanone. Must be numeric entries.
  ## Design matrix for parent phenotypes
  covM.dat <- qtlnet:::pull.pheno.null(cross, pheno.parents)
  ## Design matrix for additive covariates.
  addcovM.dat <- covM.dat
  tmp <- unique(c(addcov, intcov))
  if(!is.null(tmp)) {
    tmp <- qtlnet:::create.cov.matrix(cross, cov.names = tmp)
    if(length(tmp))
      addcovM.dat <- cbind(addcovM.dat, tmp)
  }
  ## Design matrix for interactive covariates.
  intcovM.dat <- qtlnet:::create.cov.matrix(cross,cov.names=intcov)

  ## Phenotype matrices of actual data. Can be mix of numeric and factor.
  addcov.dat <- qtlnet:::pull.pheno.null(cross, unique(c(addcov, intcov)))
  intcov.dat <- qtlnet:::pull.pheno.null(cross, intcov)

  ## Call to scanone is the big time commitment.
  #ss <- scanone.summary(cross, pheno.col, addcovM.dat, intcovM.dat, threshold, method)
  
  bic <- rep(NA, n.pheno)
  for(i in seq(n.pheno)) {
    if(i == 1) {
      y <- cross$pheno[, pheno.col[1]]
      ## Need only do this once if no QTL.
      form <- qtlnet:::set.dat.form(y, covM.dat, addcov.dat, intcov.dat)
      dat <- form$dat
      form <- form$form
        
      ## Fit linear model.
      fm <- lm(form, dat)
    }
    else
      dat$y <- cross$pheno[, pheno.col[i]]
        
    ## Record BIC.
    bic[i] <- AIC(update(fm, data = dat), k = log(length(y)))[1]
  }
  bic
}

#################################################################
bic.qtlnet.pheno <- function(cross, pheno.col, threshold,
addcov=NULL, intcov=NULL,
max.parents = 3,
parents = parents.qtlnet(pheno.col, max.parents),
verbose = TRUE,
...)
{
    ## This is modified from bic.qtlnet in R/qtlnet.
    ## This modified function works when there is no genotype information,
    ## hence, it does not run QTL mapping.
    
    ## Pre-compute BIC for many patterns.
    
    ## Calculate genotype probabilities if not already done.
    if (!("prob" %in% names(cross$geno[[1]]))) {
        warning("First running calc.genoprob.")
        cross <- calc.genoprob(cross)
    }
    
    ## Adjust phenotypes and covariates to be numeric.
    cross <- qtlnet:::adjust.pheno(cross, pheno.col, addcov, intcov)
    pheno.col <- cross$pheno.col
    pheno.names <- cross$pheno.names
    addcov <- cross$addcov
    intcov <- cross$intcov
    cross <- cross$cross
    
    n.pheno <- length(pheno.col)
    
    ## LOD threshold by phenotype.
    if(length(threshold) == 1)
    threshold <- rep(threshold, n.pheno)
    threshold <- as.list(threshold)
    if(length(threshold) != n.pheno)
    stop("threshold must have same length as pheno.col")
    
    n.scores <- sum(n.pheno - sapply(parents, length))
    update.scores <- list(code = rep("", n.scores),
    pheno.col = rep(NA, n.scores),
    bic = rep(NA, n.scores))
    k <- 0
    
    for(parent.code in names(parents)) {
        if(verbose)
        print(parent.code)
        ## Compute BIC ahead of time.
        ## This method just marches through based on parent pattern.
        ## To do this in parallel,
        pheno.cols <- which(!(seq(n.pheno) %in% parents[[parent.code]]))
        code <- qtlnet:::find.code(pheno.cols, parents[[parent.code]])
        for(i in seq(length(pheno.cols))) {
            bic <- qtlnet:::find.bic(code[i], pheno.cols[i], update.scores, ...)
            if(is.na(bic)){
                ## Modification: calc.pheno.bic instead of calc.bic in R/qtlnet
                ## calc.pheno.bic does not do QTL mapping.
                run <- calc.pheno.bic(cross, code[i], pheno.cols[i], parents[[parent.code]],
                addcov, intcov, threshold, n.pheno,
                update.scores = update.scores, ...)
                n.run <- nrow(run)
                update.scores$code[k + seq(n.run)] <- as.character(run$code)
                update.scores$pheno.col[k + seq(n.run)] <- run$pheno.col
                update.scores$bic[k + seq(n.run)] <- run$bic
                k <- k + n.run
            }
        }
    }
    update.scores
}

