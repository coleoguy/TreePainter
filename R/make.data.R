library(diversitree)
trees <- traits <- list()
for(i in 1:100){
  trees[[i]] <- trees(pars=c(3,1),
                      type="bd",
                      max.taxa = (i+99))[[1]]
  depth <- max(branching.times(trees[[i]]))
  trees[[i]]$edge.length <- trees[[i]]$edge.length / depth
  traits[[i]] <- sim.character(tree = trees[[1]],
                               pars = matrix(c(-.4,.4,.4,-.4), 2,2),
                               model = "mkn",x0 = 2)
}
library(castor)


fit <- treePaintR(tree = tree,
                  tip_states = trait,
                  qmat = qmat,
                  iter = 10000,
                  rate.classes = 7)
plotRateTree(tree = fit[[1]],
             rates = 7,
             scaled = F,
             cols=NULL,
             bg="lightgray", edge.width = 1)

plot(fit$lk.trace, type = "l",
     ylab = "lnLik",
     xlab = "generation")


