library(ape)
library(viridis)
library(phytools)

source("functions.R")
library(castor)
library(phytools)
library(diversitree)

rate.class <- 1:7
rates <- c(.7,.8,.9,1,1.1,1.2,1.3)
tips <- 100

#make tree
# set.seed(3)
tree <- rcoal(n = tips)

# make a simple model
qmat <- matrix(data = c(-.4,.4,.4,-.4),
               nrow = 2,
               ncol = 2)
colnames(qmat) <- rownames(qmat) <- c("1","2")

# get tip states
# trait <- sample(c(1,2), tips, replace = T)
trait <- sim.character(tree = tree,
                       pars = qmat,
                       model = "mkn",x0 = 2)
trait2 <- rep(1, 100)
trait2[85:100] <- sample(c(1,2), size=16, replace=T)
painted.tree <- treePaintR(tree = tree,
                           tip_states = trait,
                           qmat = qmat,
                           iter = 10000,
                           rate.class = rate.class,
                           rates = rates)
cols <- heat.colors(7)[7:1]
par(bg="gray")
plot(tree,
     show.tip.label = F,
     edge.color = cols[as.factor(painted.tree$tree$rates)],
     edge.width = 2)

plot(painted.tree$lk.trace, type = "l",
     ylab = "-loglikelihood",
     xlab = "Iteration")


### troubleshooting code
tip_states = trait
iter = 10000
