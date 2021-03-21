library(diversitree)
qmat <- matrix(c(-.4,.4,.4,-.4), 2,2)
trees <- traits <- list()
for(i in 1:100){
  trees[[i]] <- trees(pars=c(3,1),
                      type="bd",
                      max.taxa = (i+99))[[1]]
  depth <- max(branching.times(trees[[i]]))
  trees[[i]]$edge.length <- trees[[i]]$edge.length / depth
  traits[[i]] <- sim.character(tree = trees[[i]],
                               pars = qmat,
                               model = "mkn",x0 = 2)
}
library(castor)


fit <- treePaintR(tree = trees[[100]],
                  tip_states = traits[[100]],
                  qmat = qmat,
                  iter = 12000,
                  rate.classes = 9,
                  step = 9)
plotRateTree(tree = fit[[1]],
             rates = 9,
             scaled = F,
             cols=NULL,
             bg="lightgray", edge.width = 1)

plot(fit$lk.trace, type = "l",
     ylab = "lnLik",
     xlab = "generation")


