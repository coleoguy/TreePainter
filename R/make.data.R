library(diversitree)
source("functions.R")
qmat <- matrix(c(-.1,.1,.1,-.1), 2,2)
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
  sim.tree$edge.length[branches] <- sim.tree$edge.length[branches] * 2
  traits[[i]] <- as.vector(sim.character(tree = sim.tree,
                               pars = qmat,
                               model = "mkn",x0 = 2))
}
library(castor)
fit3 <- treePaintR(tree = trees[[100]],
                  tip_states = traits[[100]],
                  qmat = qmat,
                  iter = 10000,
                  rate.classes = 3,
                  step = 5)
plot(tree = fit3[[1]],
     rates = fit3$num.rates,
     scaled = F,
     cols=NULL,
     bg="lightgray", edge.width = 1)
plot(tree = fit2[[1]],
     rates = 13,
     scaled = T,
     cols=NULL,
     bg="lightgray", edge.width = 1)
plot(tree = fit3[[1]],
     rates = 13,
     scaled = T,
     cols=NULL,
     bg="lightgray", edge.width = 1)

plot(fit3$lk.trace, type = "l",
     ylab = "lnLik",
     xlab = "generation")
lines(fit$lk.trace, type = "l",col="red")
lines(fit3$lk.trace, type = "l",col="red")


