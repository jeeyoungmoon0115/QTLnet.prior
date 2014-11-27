###########################################################
## 2014/07/17
##  Copyright (C) 2014 Jee Young Moon and Brian S. Yandell
## This code is under the GNU.
##
## Functions: sim.example1.cross, geno.effect, sim.knowledge
############################################################
sim.example1.cross<-function(seed=1234){
  set.seed(seed)
  simmap <- sim.map(len=rep(100,5), n.mar=10, include.x=FALSE, eq.spacing=FALSE)
  model= rbind(c(1,50,0,0), c(2, 50,0,0), c(4,50,0,0), c(5,50,0,0))
  additive <- runif(4, min=0, max=0.5)
  dominance <- runif(4, min=0, max=0.25)
  #additive <- runif(4, min=0.2, max=0.5)
  #dominance <- runif(4, min=0.2, max=0.25)
  partialbeta <- runif(7, min=-0.5, max=0.5)
  #partialbeta <- runif(7, min=0.3, max=0.7)
  #partialbeta <- partialbeta * sample(c(1,-1), 7, replace=TRUE)
  n.ind=500
  simcross<-sim.cross(map=simmap, model=model, n.ind=n.ind, type='f2', keep.qtlgeno=TRUE)

  Y1 <- geno.effect(simcross$qtlgeno[,1], additive[1], dominance[1]) + rnorm(n.ind, 0, 1)
  Y2 <- geno.effect(simcross$qtlgeno[,2], additive[2], dominance[2]) + Y1*partialbeta[1] + rnorm(n.ind, 0,1)
  Y3 <- Y1*partialbeta[2] + rnorm(n.ind, 0,1)
  Y4 <- geno.effect(simcross$qtlgeno[,3], additive[3], dominance[3]) + Y1*partialbeta[3] + Y3 * partialbeta[4] + rnorm(n.ind, 0,1)
  Y5 <- geno.effect(simcross$qtlgeno[,4], additive[4], dominance[4]) + Y2 * partialbeta[5] + Y3*partialbeta[6] + Y4*partialbeta[7] + rnorm(n.ind, 0,1)

  simcross$pheno<-data.frame(cbind(Y1,Y2,Y3,Y4,Y5))
  simcross
}


geno.effect<-function(genotype, additive, dominance){
  ## by Cockerham's model  
  geno.add <- genotype-2
  geno.dom <- (1+geno.add)*(1-geno.add)-1/2
  return(geno.add*additive + geno.dom*dominance)
}




sim.knowledge<-function(B.N=1, delta=0.5, A){
    ## Prior knowledge matrix B: 5 x 5 matrix
    ## N_{+,-} (mu +- delta, sigma) : [0,1]-truncated normal distribution
    ## mu=0.5, sigma=0.1
    ## delta = -0.5, -0.25, 0, 0.25, 0.5
    mu=0.5
    sigma=0.1
    #delta=c(-0.5, -0.25, 0, 0.25, 0.5)
    
    if (B.N==0) return(NULL)
    if (length(delta)!=B.N) delta<-rep(delta[1], B.N)
    
    ## [0,1]-truncated normal distribution
    #require(msm)
    B <- vector("list", B.N)
    for (i in seq(B.N)){
        B[[i]]<-apply(A, c(1,2), function(x, mu, delta.i, sigma) {
            if (x==1) rtnorm(1, mean=mu+delta.i, sd=sigma, lower=0, upper=1)
            else rtnorm(1, mean=mu-delta.i, sd=sigma, lower=0, upper=1)
        }, mu, delta[i], sigma)
        
        diag(B[[i]]) <- 0
    }
    return(B)
}
