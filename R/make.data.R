library(diversitree)
source("functions.R")
qmat <- matrix(c(-.3,.3,.3,-.3), 2,2)
trees <- traits <- list()
fast.branches <- list()
for(i in 1:100){
  trees[[i]] <- trees(pars=c(3,1),
                      type="bd",
                      max.taxa = (i+99))[[1]]
  # scale tree to unit length
  depth <- max(branching.times(trees[[i]]))
  trees[[i]]$edge.length <- trees[[i]]$edge.length / depth
  working <- T
  while(working){
    hit <- sample(trees[[i]]$edge[,1], 1)
    tips <- get.descendants(node = hit,
                            tree = trees[[i]], tips.only = T)
    branches <- get.descendants(node = hit,
                            tree = trees[[i]], edge.index = T)
    if(length(tips) > .2*length(trees[[i]]$tip.label) &
       length(tips) < .8*length(trees[[i]]$tip.label)){
      working <- F
      fast.branches[[i]] <- branches
    }
  }
  sim.tree <- trees[[i]]
  sim.tree$edge.length[branches] <- sim.tree$edge.length[branches] * 5
  traits[[i]] <- as.vector(sim.character(tree = sim.tree,
                               pars = qmat,
                               model = "mkn",x0 = 2))
}
library(castor)
fit <- treePaintR(tree = trees[[100]],
                  tip_states = traits[[100]],
                  qmat = qmat,
                  iter = 10000,
                  rate.classes = 13,
                  step = 1)
par(mfcol = c(1,2))
plot(tree = fit[[1]],
     rates = fit$num.rates,
     scaled = F,
     cols=NULL,
     bg="lightgray", edge.width = 1)
plot(sim.tree, show.tip.label = F)
tiplabels(pch = 16, col = c("red", "blue")[traits[[100]]], cex = .5, offset = .01)


# plot(tree = fit2[[1]],
#      rates = 13,
#      scaled = T,
#      cols=NULL,
#      bg="lightgray", edge.width = 1)
# plot(tree = fit3[[1]],
#      rates = 13,
#      scaled = T,
#      cols=NULL,
#      bg="lightgray", edge.width = 1)
# 
# plot(fit3$lk.trace, type = "l",
#      ylab = "lnLik",
#      xlab = "generation")
# lines(fit$lk.trace, type = "l",col="red")
# lines(fit3$lk.trace, type = "l",col="red")


