library(diversitree)
qmat <- matrix(c(-.4,.4,.4,-.4), 2,2)
trees <- traits <- list()
for(i in 1:100){
  trees[[i]] <- trees(pars=c(3,1),
                      type="bd",
                      max.taxa = (i+99))[[1]]
  depth <- max(branching.times(trees[[i]]))
  trees[[i]]$edge.length <- trees[[i]]$edge.length / depth
  #TODO add a step of creating variation in rates so far every
  #thing we are doing is a null result
  traits[[i]] <- as.vector(sim.character(tree = trees[[i]],
                               pars = qmat,
                               model = "mkn",x0 = 2))
}
library(castor)
tips <- get.descendants(179, trees[[1]])
traits[[1]][tips] <- sample(1:2, length(tips), replace=T)
fit3 <- treePaintR(tree = trees[[1]],
                  tip_states = traits[[1]],
                  qmat = qmat,
                  iter = 10000,
                  rate.classes = 13,
                  step = 9)
plot(tree = fit3[[1]],
     rates = 13,
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


